"""
Display utilities for MolForge CLI.

Provides ASCII art banners, colored output, and formatting functions.
"""

import sys
from typing import Optional


# ANSI color codes
class Colors:
    """ANSI color code constants."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    
    @staticmethod
    def is_supported() -> bool:
        """Check if terminal supports colors."""
        return hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()


def colorize(text: str, color: str) -> str:
    """
    Colorize text if terminal supports it.
    
    Args:
        text: Text to colorize
        color: Color code from Colors class
        
    Returns:
        Colorized text if supported, plain text otherwise
    """
    if Colors.is_supported():
        return f"{color}{text}{Colors.RESET}"
    return text


def print_banner():
    """Print MolForge ASCII art banner (full box; kept for optional use)."""
    from .. import __version__
    _w = 63
    tagline = "Molecular data processing, curation & conformer generation".center(_w)
    version = f"v{__version__}".center(_w)
    banner = f"""
╔═══════════════════════════════════════════════════════════════╗
║   ███╗   ███╗ ██████╗ ██╗     ███████╗ ██████╗ ██████╗  ██████╗ ███████╗   ║
║   ████╗ ████║██╔═══██╗██║     ██╔════╝██╔═══██╗██╔══██╗██╔════╝ ██╔════╝   ║
║   ██╔████╔██║██║   ██║██║     █████╗  ██║   ██║██████╔╝██║  ███╗█████╗     ║
║   ██║╚██╔╝██║██║   ██║██║     ██╔══╝  ██║   ██║██╔══██╗██║   ██║██╔══╝     ║
║   ██║ ╚═╝ ██║╚██████╔╝███████╗██║     ╚██████╔╝██║  ██║╚██████╔╝███████╗   ║
║   ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝      ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚══════╝   ║
║{tagline}║
║{version}║
╚═══════════════════════════════════════════════════════════════╝
"""
    print(colorize(banner, Colors.CYAN))


def print_simple_banner():
    """Print the simplified MolForge banner."""
    from .. import __version__
    art = r"""
    __  ___      __   ______
   /  |/  /___  / /  / ____/___  _________ ____
  / /|_/ / __ \/ /  / /_  / __ \/ ___/ __ `/ _ \
 / /  / / /_/ / /  / __/ / /_/ / /  / /_/ /  __/
/_/  /_/\____/_/  /_/    \____/_/   \__, /\___/
                                   /____/
"""
    tagline = f"    Molecular data processing, curation & conformer generation  v{__version__}"
    print(colorize(art + tagline, Colors.CYAN))


def print_diagram():
    """Print an architecture overview built dynamically from the registry."""
    from ..configuration.registry import ActorRegistry
    from ..backends.registry import BackendRegistry

    lines = [
        "┌─────────────────────────────────────────────┐",
        "│           MolForge Architecture              │",
        "└─────────────────────────────────────────────┘",
        "",
        "Core steps (executed in the order configured):",
    ]
    for step in ActorRegistry.get_core_steps():
        info = ActorRegistry.get_actor_info(step)
        actor_cls = info['class'].__name__ if info else '?'
        backends = BackendRegistry.list_backends(step).get(step, [])
        suffix = f"  [backends: {', '.join(sorted(backends))}]" if backends else ""
        lines.append(f"  {step:<14} → {actor_cls}{suffix}")

    plugins = ActorRegistry.get_plugin_steps()
    if plugins:
        lines += ["", f"Plugins: {', '.join(plugins)}"]

    lines += ["", "Run 'molforge info <step>' for a step's parameters."]
    print("\n".join(lines))


def print_success(message: str):
    """Print success message in green with checkmark."""
    print(f"{colorize('✓', Colors.GREEN)} {message}")


def print_error(message: str):
    """Print error message in red with X."""
    print(f"{colorize('✗', Colors.RED)} {message}", file=sys.stderr)


def print_warning(message: str):
    """Print warning message in yellow with warning symbol."""
    print(f"{colorize('⚠', Colors.YELLOW)} {message}")


def print_info(message: str):
    """Print info message in blue."""
    print(f"{colorize('ℹ', Colors.BLUE)} {message}")


def print_step(step: int, total: int, name: str, rows: int, time: float):
    """
    Print pipeline step completion.
    
    Args:
        step: Current step number
        total: Total number of steps
        name: Step name
        rows: Number of rows processed
        time: Time taken in seconds
    """
    status = colorize('✓', Colors.GREEN)
    time_str = format_time(time)
    print(f"  {status} [{step}/{total}] {name:<16} | {rows:>5} rows | {time_str}")


def format_time(seconds: float) -> str:
    """
    Format time duration in human-readable format.
    
    Args:
        seconds: Duration in seconds
        
    Returns:
        Formatted time string (e.g., "1.2s", "2m 30s", "1h 15m")
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes, secs = divmod(seconds, 60)
        return f"{int(minutes)}m {secs:.0f}s"
    else:
        hours, remainder = divmod(seconds, 3600)
        minutes, secs = divmod(remainder, 60)
        return f"{int(hours)}h {int(minutes)}m {secs:.0f}s"


def progress_bar(current: int, total: int, width: int = 40) -> str:
    """
    Generate a simple progress bar.
    
    Args:
        current: Current progress value
        total: Total value
        width: Width of progress bar in characters
        
    Returns:
        Progress bar string
    """
    filled = int(width * current / total) if total > 0 else 0
    bar = '━' * filled + '─' * (width - filled)
    percent = 100 * current / total if total > 0 else 0
    return f"{bar} {percent:.0f}% | {current}/{total}"


def print_section(title: str):
    """Print section header."""
    print(f"\n{colorize(title, Colors.BOLD)}")
    print("─" * len(title))


def print_table_row(col1: str, col2: str, width1: int = 20):
    """Print a simple two-column table row."""
    print(f"  {col1:<{width1}} {col2}")
