
# =============================================================================
# Logging Setup (logger.py)
# =============================================================================

"""
Enhanced logging system with Rich console integration.
Provides structured logging with progress bars and colored output.
"""

import logging
import sys
from pathlib import Path

try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.progress import Progress, TaskID, BarColumn, TextColumn, TimeElapsedColumn
    from rich.table import Table
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

class GenomeAnalysisLogger:
    """Enhanced logger with Rich integration for better console output."""
    
    def __init__(self, level=logging.INFO, log_file=None, use_rich=True):
        self.use_rich = use_rich and RICH_AVAILABLE
        self.console = Console() if self.use_rich else None
        self.progress = None
        self.current_tasks = {}
        
        # Setup logging
        self.logger = logging.getLogger('genome_analysis')
        self.logger.setLevel(level)
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Setup handlers
        if self.use_rich:
            handler = RichHandler(console=self.console, show_time=True, show_path=False)
            handler.setFormatter(logging.Formatter('%(message)s'))
        else:
            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            ))
        
        self.logger.addHandler(handler)
        
        # File handler if specified
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            ))
            self.logger.addHandler(file_handler)
    
    def info(self, message):
        """Log info message."""
        self.logger.info(message)
    
    def warning(self, message):
        """Log warning message."""
        self.logger.warning(message)
    
    def error(self, message):
        """Log error message."""
        self.logger.error(message)
    
    def debug(self, message):
        """Log debug message."""
        self.logger.debug(message)
    
    def section_header(self, title, width=80):
        """Print a formatted section header."""
        if self.use_rich:
            self.console.print(f"\n{'=' * width}")
            self.console.print(f"{title}", style="bold blue")
            self.console.print(f"{'=' * width}")
        else:
            print(f"\n{'=' * width}")
            print(f"{title}")
            print(f"{'=' * width}")
    
    def create_progress_bar(self, description="Processing"):
        """Create a new progress bar."""
        if self.use_rich:
            self.progress = Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TextColumn("({task.completed}/{task.total})"),
                TimeElapsedColumn(),
                console=self.console
            )
            self.progress.start()
            return self.progress.add_task(description, total=100)
        else:
            return None
    
    def update_progress(self, task_id, advance=1, description=None):
        """Update progress bar."""
        if self.progress and task_id is not None:
            if description:
                self.progress.update(task_id, description=description, advance=advance)
            else:
                self.progress.advance(task_id, advance)
    
    def finish_progress(self):
        """Finish and cleanup progress bar."""
        if self.progress:
            self.progress.stop()
            self.progress = None
    
    def create_summary_table(self, title, data):
        """Create a formatted summary table."""
        if self.use_rich:
            table = Table(title=title, show_header=True, header_style="bold magenta")
            table.add_column("Metric", style="dim")
            table.add_column("Value", justify="right")
            
            for key, value in data.items():
                if isinstance(value, float):
                    table.add_row(key, f"{value:.3f}")
                else:
                    table.add_row(key, str(value))
            
            self.console.print(table)
        else:
            print(f"\n{title}")
            print("-" * len(title))
            for key, value in data.items():
                if isinstance(value, float):
                    print(f"{key}: {value:.3f}")
                else:
                    print(f"{key}: {value}")

# Global logger instance
_logger = None

def setup_logger(level=logging.INFO, log_file=None, use_rich=True):
    """Setup the global logger instance."""
    global _logger
    _logger = GenomeAnalysisLogger(level=level, log_file=log_file, use_rich=use_rich)
    return _logger

def get_logger():
    """Get the global logger instance."""
    global _logger
    if _logger is None:
        _logger = setup_logger()
    return _logger
