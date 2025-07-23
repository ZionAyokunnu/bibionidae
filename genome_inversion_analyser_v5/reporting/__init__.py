# =============================================================================
# Reporting Package (__init__.py)
# =============================================================================

"""
Comprehensive reporting system for genome analysis results.
Provides formatted reports and statistical summaries.
"""

from .report_generator import ReportGenerator
from .statistics_formatter import StatisticsFormatter

__all__ = [
    'ReportGenerator',
    'StatisticsFormatter'
]