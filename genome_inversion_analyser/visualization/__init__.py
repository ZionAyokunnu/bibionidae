"""Visualization module for genome inversion analyzer"""

from .plots import (
    create_enhanced_visualizations,
    add_synteny_block_lines,
    create_busco_synteny_dotplot,
    create_ortholog_quality_plots,
    create_synteny_block_plots,
    create_inversion_landscape_plot,
    create_rearrangement_network_plot,
    create_analysis_summary_dashboard,
    create_statistics_summary_plot,
    create_quality_summary_plot,
    create_chromosome_mapping_overview,
    create_synteny_summary_plot,
    create_inversion_summary_plot,
    create_rearrangement_summary_plot,
    create_method_summary_plot
)
from .publication_plots import (
    PublicationPlotGenerator,
    create_publication_plots,
    create_curved_synteny_plot,
    create_annotated_phylogeny
)
from .syri_integration import SyRIIntegrator

__all__ = [
    'create_enhanced_visualizations',
    'add_synteny_block_lines',
    'create_busco_synteny_dotplot',
    'create_ortholog_quality_plots',
    'create_synteny_block_plots',
    'create_inversion_landscape_plot',
    'create_rearrangement_network_plot',
    'create_analysis_summary_dashboard',
    'create_statistics_summary_plot',
    'create_quality_summary_plot',
    'create_chromosome_mapping_overview',
    'create_synteny_summary_plot',
    'create_inversion_summary_plot',
    'create_rearrangement_summary_plot',
    'create_method_summary_plot'
]