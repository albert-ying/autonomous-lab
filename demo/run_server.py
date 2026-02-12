#!/usr/bin/env python3
"""Start the Autonomous Lab web UI pointing at the mock project."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

DEMO_DIR = "/tmp/autolab_demo_full"

def main():
    from autonomous_lab.web.main import get_web_ui_manager
    import time

    manager = get_web_ui_manager()
    manager.lab_project_dir = DEMO_DIR

    if manager.server_thread is None or not manager.server_thread.is_alive():
        manager.start_server()
        time.sleep(1.5)

    url = f"{manager.get_server_url()}/lab"
    print(f"\n  Lab UI: {url}\n")
    manager.open_browser(url)

    # Keep running
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\nShutting down...")

if __name__ == "__main__":
    main()
