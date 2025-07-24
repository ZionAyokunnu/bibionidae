"""Pipeline module for genome inversion analyzer"""

from .multi_species import MultiSpeciesPipeline, create_multi_species_config

__all__ = [
    'MultiSpeciesPipeline',
    'create_multi_species_config'
]