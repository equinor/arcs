import os
import pathlib
import time
from multiprocessing import Process, Condition
import setproctitle
import webview
import warnings 
from arcs.dash_app.domino import terminate_when_process_dies
from arcs.dash_app.server import start_dash


def get_directory_size(directory: pathlib.Path) -> int:
    """Returns the total size of the directory in bytes."""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            file_path = pathlib.Path(dirpath) / filename
            total_size += file_path.stat().st_size
    return total_size


def dataset_location() -> str:
    """Return location string of the directory with dataset to be used."""
    small_dataset = os.path.join(os.path.dirname(__file__),'data/')
    large_dataset = os.path.join(os.path.dirname(__file__),'data/large_dataset/')

    # check if large dataset exists and that git-lfs pulled has happened
    directory = pathlib.Path(large_dataset)
    if directory.exists() and directory.is_dir():
        MANY_MEGABYTES = 100*1024*1024
        if get_directory_size(large_dataset) > MANY_MEGABYTES:
            return large_dataset
    
    # otherwise use the small data set
    return small_dataset


def start():
    warnings.simplefilter('ignore')
    port = int(os.getenv("PORT", "8050"))
    host = os.getenv("HOST", "127.0.0.1")
    this_dir, this_filename = os.path.split(__file__)
    file_location = dataset_location()
    server_is_started = Condition()

    # Set the process title.
    setproctitle.setproctitle('arcs-0.1.0')

    # Spawn the dash process.
    p = Process(target=start_dash, args=(host, port, server_is_started, file_location))
    p.start()
    # If the dash process dies, follow along.
    terminate_when_process_dies(p)

    # Wait until dash process is ready.
    with server_is_started:
        server_is_started.wait()

    time.sleep(0.2)

    # Create the webview.
    webview.create_window('ARCS 1.4.0', f'http://{host}:{port}',
                          width=1000, 
                          height=1000)
    webview.start()

    # Reached when window is closed.
    p.terminate()
    exit(0)

if __name__ == '__main__':
    start()
