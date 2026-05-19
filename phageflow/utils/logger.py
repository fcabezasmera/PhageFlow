"""PhageFlow logging utilities using Rich."""

from datetime import datetime
from pathlib import Path
from rich.console import Console
from rich.theme import Theme
from rich.panel import Panel
from rich.text import Text

_THEME = Theme({
    "info":    "cyan",
    "ok":      "bold green",
    "warn":    "bold yellow",
    "error":   "bold red",
    "step":    "bold cyan",
    "muted":   "dim white",
})

console = Console(theme=_THEME)


def log_info(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    console.print(f"[muted]\\[{ts}][/muted] [info]INFO [/info]  {msg}")


def log_ok(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    console.print(f"[muted]\\[{ts}][/muted] [ok]OK   [/ok]   {msg}")


def log_warn(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    console.print(f"[muted]\\[{ts}][/muted] [warn]WARN [/warn]  {msg}")


def log_error(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    console.print(f"[muted]\\[{ts}][/muted] [error]ERROR[/error] {msg}")


def log_step(title: str) -> None:
    console.rule(f"[step]{title}[/step]", style="cyan")


def log_header(version: str) -> None:
    text = Text()
    text.append("PhageFlow", style="bold cyan")
    text.append(f" v{version}", style="cyan")
    text.append("\nModular bacteriophage genomics pipeline", style="dim white")
    console.print(Panel(text, border_style="cyan", padding=(0, 2)))


def setup_file_log(log_path: Path) -> None:
    """Write a duplicate log to a file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
