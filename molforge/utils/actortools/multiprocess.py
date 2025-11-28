from typing import Tuple, List, Optional, Any, Callable

import gc

import multiprocessing as mp
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

from ...actors.params.base import BaseParams
from ...configuration.logger import PipelineLogger
from ..constants import MAX_CHUNK_SIZE


def calculate_chunk_params(data_length: int, max_chunk_size: int = MAX_CHUNK_SIZE) -> Tuple[int, int, int]:
    """
    Calculate optimal chunking parameters for multiprocessing.

    Args:
        data_length: Total number of items to process
        max_chunk_size: Maximum items per chunk (default from constants)

    Returns:
        Tuple of (n_processes, chunk_size, n_chunks)
        - n_processes: Number of worker processes (CPU count - 1, minimum 1)
        - chunk_size: Items per chunk (evenly distributed, capped at max_chunk_size)
        - n_chunks: Total number of chunks needed
    """
    n_processes = max(1, mp.cpu_count() - 1)
    chunk_size = min((data_length + n_processes - 1) // n_processes, max_chunk_size)
    n_chunks = (data_length + chunk_size - 1) // chunk_size
    return n_processes, chunk_size, n_chunks


def multiprocess_worker(chunk: List[Any],
                       actor_class: type,
                       actor_method: str,  # Method name as string
                       actor_params: BaseParams,
                       logger: Optional[PipelineLogger] = None,
                       method_kwargs: dict = None) -> List[Tuple]:
    """
    Generalized worker function for multiprocessing any actor class method.

    Args:
        chunk (List[Any]): List of individual serializable inputs to process.
        actor_class (type): The actor class to instantiate for this worker.
        actor_method (str): Name of the method to call on the actor instance.
        actor_params (BaseParams): Instance of parameters for this actor class.
        logger (Optional[PipelineLogger], optional): Logger instance for logging. Defaults to None.
        method_kwargs (dict, optional): Keyword arguments to pass to the method. Defaults to None.

    Returns:
        List[Tuple]: chunk results.
    """    
    import traceback
    import sys
    import gc
    import os
    import time
    import signal
    
    if method_kwargs is None:
        method_kwargs = {}
    
    worker_pid = os.getpid()
    start_time = time.time()
    processed_count = 0
    
    
    # Enhanced error tracking
    def log_error(message, level='ERROR'):
        try:
            if logger:
                logger.log(level.lower(), '[MULTIP]', message)

            print(f"WORKER {worker_pid}: {message}")
        except:
            print(f"WORKER {worker_pid}: {message}")
    
    # Signal handlers
    def signal_handler(signum, frame):
        signal_names = {
            signal.SIGTERM: 'SIGTERM',
            signal.SIGINT: 'SIGINT',
            signal.SIGSEGV: 'SIGSEGV', 
            signal.SIGABRT: 'SIGABRT'
        }
        signal_name = signal_names.get(signum, f'Signal {signum}')
        error_msg = f"Received {signal_name} after {time.time() - start_time:.2f}s, processed {processed_count}/{len(chunk)} items"
        log_error(f"SIGNAL: {error_msg}")
        
        try:
            print(f"\nStack trace when {signal_name} received:", file=sys.stderr)
            traceback.print_stack(frame, file=sys.stderr)
        except:
            pass
            
        raise KeyboardInterrupt(error_msg)
    
    # Register signal handlers
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    if hasattr(signal, 'SIGSEGV'):
        signal.signal(signal.SIGSEGV, signal_handler)
    if hasattr(signal, 'SIGABRT'):
        signal.signal(signal.SIGABRT, signal_handler)
    
    try:
        # Force garbage collection at start
        gc.collect()
        
        log_error(f"Starting with {len(chunk)} items using {actor_class.__name__}.{actor_method}", level='debug')

        # Initialize actor with logger (suppress init logs for workers)
        try:
            actor = actor_class(actor_params, logger=logger, suppress_init_logs=True)
        except Exception as e:
            error_msg = f"Failed to initialize {actor_class.__name__}: {type(e).__name__}: {str(e)}\nTraceback:\n{traceback.format_exc()}"
            log_error(f"INIT ERROR: {error_msg}")
            return [(None, False, f'init_error: {error_msg}') for _ in chunk]
        
        # Get the method to call
        try:
            method = getattr(actor, actor_method)
            if not callable(method):
                raise AttributeError(f"{actor_method} is not callable")
        except AttributeError as e:
            error_msg = f"Method {actor_method} not found on {actor_class.__name__}: {str(e)}"
            log_error(f"METHOD ERROR: {error_msg}")
            return [(None, False, f'method_error: {error_msg}') for _ in chunk]
        
        results = []
        
        for i, item in enumerate(chunk):
            processed_count = i
            
            try:
                # Call the method with the item and any additional kwargs
                result = method(item, **method_kwargs)
                results.append(result)
                
            except Exception as e:
                # Detailed error for this specific item
                error_details = f"Item {i} failed: {type(e).__name__}: {str(e)}"
                results.append((None, False, f'item_error: {error_details}'))
                
                # Log full traceback for debugging
                log_error(f"Item {i} error:\n{traceback.format_exc()}", level='debug')
            
            # Periodic garbage collection for large chunks
            if i % 1000 == 0 and i > 0:
                gc.collect()
                # log_error(f"Processed {i}/{len(chunk)} items", level='debug') # too verbose

        processed_count = len(chunk)
        
        # Final cleanup and validation
        gc.collect()
        
        if len(results) != len(chunk):
            missing = len(chunk) - len(results)
            error_msg = f"Result count mismatch: expected {len(chunk)}, got {len(results)}"
            log_error(f"COUNT ERROR: {error_msg}")
            results.extend([(None, False, f'missing_result: {error_msg}') for _ in range(missing)])
        
        runtime = time.time() - start_time
        log_error(f"Completed successfully: {len(results)} results in {runtime:.2f}s", level='debug')
        
        return results
        
    except MemoryError as e:
        error_msg = f"Memory error after {time.time() - start_time:.2f}s: {str(e)} (processed {processed_count}/{len(chunk)})"
        log_error(f"MEMORY ERROR: {error_msg}")
        return [(None, False, f'memory_error: {error_msg}') for _ in chunk]
        
    except KeyboardInterrupt as e:
        error_msg = f"Interrupted after {time.time() - start_time:.2f}s: {str(e)} (processed {processed_count}/{len(chunk)})"
        log_error(f"INTERRUPTED: {error_msg}")
        return [(None, False, f'interrupted: {error_msg}') for _ in chunk]
        
    except Exception as e:
        tb_str = traceback.format_exc()
        runtime = time.time() - start_time
        error_msg = f"Unexpected error after {runtime:.2f}s: {type(e).__name__}: {str(e)} (processed {processed_count}/{len(chunk)})\nFull traceback:\n{tb_str}"
        log_error(f"UNEXPECTED ERROR: {error_msg}")
        return [(None, False, f'worker_error: {error_msg}') for _ in chunk]
    
    finally:
        try:
            gc.collect()
        except:
            pass