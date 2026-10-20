"""Launch Venison as a desktop app.

Starts a local Panel server and shows it in a native window (pywebview), or in
your default browser if pywebview isn't available or --no-webview is given.

    python launch_venison.py [--port 5007] [--no-webview] [--no-browser]
"""

import argparse
import logging
import multiprocessing as mp
import os
import socket
import sys
import threading
import time
import warnings
import webbrowser

LOADING_HTML = """<!DOCTYPE html><html><head><style>
body{margin:0;height:100vh;display:flex;align-items:center;justify-content:center;
background:#FAF9F7;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.w{font-size:34px;font-weight:700;letter-spacing:-.02em;color:#26211E}
.s{margin-top:10px;font-size:13px;color:#7A706A;text-align:center}
</style></head><body><div><div class="w">Venison</div>
<div class="s">Starting…</div></div></body></html>"""


def wait_for_server(port: int, timeout: float = 30) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        with socket.socket() as sock:
            sock.settimeout(1)
            if sock.connect_ex(("localhost", port)) == 0:
                return True
        time.sleep(0.25)
    return False


def port_in_use(port: int) -> bool:
    with socket.socket() as sock:
        return sock.connect_ex(("localhost", port)) == 0


class _BenignBokehNoise(logging.Filter):
    """Drop warnings Bokeh emits when a tab's plots are replaced mid-update."""

    def filter(self, record):
        return "Dropping a patch because it contains a previously known reference" \
            not in record.getMessage()


def serve(port: int) -> None:
    import panel as pn

    warnings.filterwarnings("ignore", message="reference already known")
    logging.getLogger().addFilter(_BenignBokehNoise())
    for handler in logging.getLogger().handlers:
        handler.addFilter(_BenignBokehNoise())

    from venison.app import configure_panel, create_app

    # Serve JS/CSS from this process instead of a CDN, so the app works offline.
    os.environ.setdefault("BOKEH_RESOURCES", "inline")
    pn.config.resources = "server"
    configure_panel()
    pn.serve(create_app, port=port, address="localhost", show=False, title="Venison",
             websocket_origin=f"localhost:{port}", threaded=False, log_level="error",
             session_token_expiration=3600)


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the Venison DEER analysis app")
    parser.add_argument("--port", type=int, default=5007)
    parser.add_argument("--no-webview", action="store_true",
                        help="open in your default browser instead of a native window")
    parser.add_argument("--no-browser", action="store_true",
                        help="only start the server and print its URL")
    args = parser.parse_args()

    # The analysis runs multiprocessing pools; spawn keeps workers from
    # inheriting the server's threads and sockets on every platform.
    mp.freeze_support()
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass

    if port_in_use(args.port):
        sys.exit(f"Port {args.port} is already in use; try --port {args.port + 1}")

    from venison.core import log_capture
    log_capture.install()

    url = f"http://localhost:{args.port}"
    threading.Thread(target=serve, args=(args.port,), daemon=True).start()

    webview = None
    if not (args.no_webview or args.no_browser):
        try:
            import webview
        except Exception:
            print("pywebview unavailable; opening in your browser instead.")

    if webview is None:
        if not wait_for_server(args.port):
            sys.exit("Server failed to start")
        if args.no_browser:
            print(f"Venison is running at {url}")
        else:
            webbrowser.open(url)
        print("Press Ctrl+C to stop.")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            os._exit(0)

    window = webview.create_window("Venison", html=LOADING_HTML, width=1440, height=900,
                                   min_size=(1000, 640))

    def load():
        if wait_for_server(args.port):
            window.load_url(url)
        else:
            window.load_html("<h2 style='font-family:sans-serif;text-align:center;"
                             "margin-top:40vh'>Venison failed to start</h2>")

    threading.Thread(target=load, daemon=True).start()
    webview.start()
    os._exit(0)


if __name__ == "__main__":
    main()
