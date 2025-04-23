from pathlib import Path
from typing import Optional
import logging
import signal
import multiprocessing as mp
import os

def format_short_path(path_obj: Path, 
                      max_levels: int = 3) -> str:
    """
    Formats a Path object to show only the last 'max_levels'
    components, prepended with '...' if the path was longer.
    This is useful for making logging messages more readable.

    Parameters:
    -----------
        path_obj: Path object to format.
        max_levels: The maximum number of path components (levels) to display.
    Returns:
    --------
        A string representation of the potentially shortened path.
    """
    # Ensure input is a Path object
    if not isinstance(path_obj, Path):
        try:
            path_obj = Path(path_obj)
        except TypeError:
            # Handle cases where input cannot be converted to Path
            return str(path_obj)
    
    parts = path_obj.parts
    if len(parts) <= max_levels:
        # If the path has max_levels or fewer components, return it as is
        return str(path_obj)
    else:
        # If the path is longer, take the last 'max_levels' components
        relevant_parts = parts[-max_levels:]
        # Join these parts with the OS-specific separator ('/' or '\')
        # Prepend with '...' to indicate truncation
        # Using os.sep ensures cross-platform compatibility
        return "..." + os.sep + os.path.join(*relevant_parts)
    
# ====================================================
# Functions to setup logging and interruption handling 
# ====================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('R2_dmplex')

_interrupt_event = None

def setup_interrupt_handling():
    """Set up interrupt handling for the process."""
    global _interrupt_event
    if _interrupt_event is None:
        _interrupt_event = mp.Event()
    
    def signal_handler(sig):
        logger.warning(f"Received signal {sig}, initiating graceful shutdown")
        if _interrupt_event:
            _interrupt_event.set()
    
    # Register signal handlers
    if mp.current_process().name == 'MainProcess':
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
    else:
        # Worker processes should use a different approach
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)


def setup_logging(
        prefix: str, 
        log_dir: Optional[Path] = None) -> None:
    """
    Set up logging configuration based on command line arguments.
    
    Parameters:
    -----------
    prefix : str
        This is used to name the log file.
    log_dir : Optional[Path]
        Directory where log files should be stored, if None use stdout
    """
    log_level = logging.INFO
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add file handler if log_dir is specified, else use console only
    if log_dir:
        # Create log directory and file if they doesn't exist
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file_path = log_dir / f"{prefix}.log"
        # If the file already exists, delete it to start fresh
        if log_file_path.is_file():
            try:
                log_file_path.unlink()
            except OSError as e:
                logger.error(f"Failed to delete existing log file: {e}")
                return
        # Create the log file
        try:
            log_file_path.touch(exist_ok=True)
        except OSError as e:
            logger.error(f"Failed to create log file: {e}")
            return
        
        # Create file handler
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
        
        # Also add console handler for immediate feedback
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)

        logger.info(f"Logging to file: {log_file_path}")
    else:
        # Add console handler only
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
