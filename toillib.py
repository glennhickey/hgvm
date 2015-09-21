"""
toillib.py: useful extras for Toil scripts.

Includes real-time-logging, input retrieval from file/S3/Azure, and output
deposit to file/S3/Azure
"""

import sys, os, os.path, json, collections, logging, logging.handlers
import SocketServer, struct, socket, threading, tarfile, shutil

def robust_makedirs(directory):
    """
    Make a directory when other nodes may be trying to do the same on a shared
    filesystem.
    
    """

    if not os.path.exists(directory):
        try:
            # Make it if it doesn't exist
            os.makedirs(directory)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass
            
    # Make sure it exists and is a directory
    assert(os.path.exists(directory) and os.path.isdir(directory))

class LoggingDatagramHandler(SocketServer.DatagramRequestHandler):
    """
    Receive logging messages from the jobs and display them on the master.
    
    Uses length-prefixed JSON message encoding.
    """
    
    def handle(self):
        """
        Handle messages coming in over self.connection.
        
        Messages are 4-byte-length-prefixed JSON-encoded logging module records.
        """
        
        while True:
            # Loop until we run out of messages
        
            # Parse the length
            length_data = self.rfile.read(4)
            if len(length_data) < 4:
                # The connection was closed, or we didn't get enough data
                # TODO: complain?
                break
                
            # Actually parse the length
            length = struct.unpack(">L", length_data)[0]
            
            # This is where we'll put the received message
            message_parts = []
            length_received = 0
            while length_received < length:
                # Keep trying to get enough data
                part = self.rfile.read(length - length_received)
                
                length_received += len(part)
                message_parts.append(part)
                
            # Stitch it all together
            message = "".join(message_parts)

            try:
            
                # Parse it as JSON
                message_attrs = json.loads(message)
                
                # Fluff it up into a proper logging record
                record = logging.makeLogRecord(message_attrs)
            except:
                logging.error("Malformed record")
                
            # TODO: do log level filtering
            logging.getLogger("remote").handle(record)
            
class JSONDatagramHandler(logging.handlers.DatagramHandler):
    """
    Send logging records over UDP serialized as JSON.
    """
    
    def makePickle(self, record):
        """
        Actually, encode the record as length-prefixed JSON instead.
        """
        
        json_string = json.dumps(record.__dict__)
        length = struct.pack(">L", len(json_string))
        
        return length + json_string
        
class RealTimeLogger(object):
    """
    All-static class for getting a logger that logs over UDP to the master.
    """
    
    # Also the logger
    logger = None
    
    # The master keeps a server and thread
    logging_server = None
    server_thread = None
  
    @classmethod
    def start_master(cls):
        """
        Start up the master server and put its details into the options
        namespace.
        
        """
        
        logging.basicConfig(level=logging.DEBUG)
    
        # Start up the logging server
        cls.logging_server = SocketServer.ThreadingUDPServer(("0.0.0.0", 0),
            LoggingDatagramHandler)
            
        # Set up a thread to do all the serving in the background and exit when we
        # do
        cls.server_thread = threading.Thread(
            target=cls.logging_server.serve_forever)
        cls.server_thread.daemon = True
        cls.server_thread.start()
        
        # Set options for logging in the class and the options namespace
        # Save them in the environment so they get sent out to jobs
        os.environ["RT_LOGGING_HOST"] = socket.getfqdn()
        os.environ["RT_LOGGING_PORT"] = str(
            cls.logging_server.server_address[1])
        
        
    @classmethod
    def stop_master(cls):
        """
        Stop the server on the master.
        
        """
        
        cls.logging_server.shutdown()
        cls.server_thread.join()
  
    @classmethod
    def get(cls):
        """
        Get the logger that logs to master.
        
        Note that if the master logs here, you will see the message twice,
        since it still goes to the normal log handlers too.
        """
        
        if cls.logger is None:
            # Only do the setup once, so we don't add a handler every time we
            # log
            cls.logger = logging.getLogger('realtime')
            cls.logger.setLevel(logging.DEBUG)
            cls.logger.addHandler(JSONDatagramHandler(
                os.environ["RT_LOGGING_HOST"], os.environ["RT_LOGGING_PORT"]))
        
        return cls.logger

def write_global_directory(file_store, path, cleanup=False):
    """
    Write the given directory into the file store, and return an ID that can be
    used to retrieve it. Writes the files in the directory and subdirectories
    into a tar file in the file store.

    Does not preserve the name or permissions of the given directory (only of
    its contents).

    If cleanup is true, directory will be deleted from the file store when this
    job and its follow-ons finish.
    
    """
    
    with file_store.writeGlobalFileStream(cleanup=cleanup) as (file_handle,
        file_id):
        # We have a stream, so start taring into it
    
        with tarfile.open(fileobj=file_handle, mode="w|gz") as tar:
            # Open it for streaming-only write (no seeking)
            
            # We can't just add the root directory, since then we wouldn't be
            # able to extract it later with an arbitrary name.
            
            for file_name in os.listdir(path):
                # Add each file in the directory to the tar
                tar.add(file_name)
                
    # Spit back the ID to use to retrieve it
    return directory_id
        
def read_global_directory(file_store, directory_id, path):
    """
    Reads a directory with the given tar file id from the global file store and
    recreates it at the given path.
    
    The given path, if it exists, must be a directory.
    
    Do not use to extract untrusted directories, since they could sneakily plant
    files anywhere on the filesystem.
    
    """
    
    # Make the path
    robust_makedirs(path)
    
    with file_store.readGlobalFileStream(directory_id) as file_handle:
        # We need to pull files out of this tar stream
    
        with tarfile.open(fileobj=file_handle, mode="r|*") as tar:
            # Open it for streaming-only read (no seeking)
            
            # We need to extract the whole thing into that new directory
            tar.extractall(path)
            
def read_input_file(input_path, local_path):
    """
    Read an input file from wherever the input comes from and send it to the
    given path.
    
    """
    
    # For now this is local filesystem. Just make a symlink
    os.symlink(os.path.abspath(input_path), local_path)
    
def list_input_directory(input_path):
    """
    Yields each of the subdirectories and files in the given input path, non-
    recursively.
    
    Gives bare file/directory names with no paths.
    
    """
    
    # Just use local files for now
    
    for item in os.listdir(input_path):
        yield item
    
def write_output_file(local_path, output_path):
    """
    Save the given local file to the given output path. No output directory
    needs to exist already.
    
    """
    
    
    parent_dir = os.path.split(output_path)[0]
    if parent_dir != "":
        # Make sure the directory it goes in exists.
        robust_makedirs(parent_dir)
    
    # These are small so we just make copies
    shutil.copy2(local_path, output_path)
    





















