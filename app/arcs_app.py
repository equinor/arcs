import os
import time
from multiprocessing import Process, Condition
import setproctitle
import webview
import warnings 
from arcs.dash_app.domino import terminate_when_process_dies
from arcs.dash_app.server_2 import start_dash

def start():
    warnings.simplefilter('ignore')
    port = int(os.getenv("PORT", "8050"))
    host = os.getenv("HOST", "127.0.0.1")
    this_dir, this_filename = os.path.split(__file__)
    file_location = os.path.join(os.path.dirname(__file__),'data/')
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
