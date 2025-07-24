"""
File Registry System for Genome Inversion Analyzer
Manages standardized outputs, dependency tracking, and format conversion
"""

import json
import pickle
import hashlib
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from datetime import datetime
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class FileRegistry:
    """
    Centralized registry for managing analysis outputs with dependency tracking
    """
    
    def __init__(self, base_dir: Union[str, Path], project_name: str = "genome_analysis"):
        self.base_dir = Path(base_dir)
        self.project_name = project_name
        self.registry_file = self.base_dir / "registry.json"
        self.registry = self._load_or_create_registry()
        
        # Create standardized directory structure
        self.directories = {
            'data': self.base_dir / 'data',
            'cache': self.base_dir / 'cache', 
            'exports': self.base_dir / 'exports',
            'bed': self.base_dir / 'exports' / 'bed',
            'gff': self.base_dir / 'exports' / 'gff',
            'json': self.base_dir / 'exports' / 'json',
            'fasta': self.base_dir / 'exports' / 'fasta',
            'plots': self.base_dir / 'plots',
            'reports': self.base_dir / 'reports'
        }
        
        self._create_directories()
        
    def _create_directories(self):
        """Create standardized directory structure"""
        for name, path in self.directories.items():
            path.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Created directory: {path}")
    
    def _load_or_create_registry(self) -> Dict:
        """Load existing registry or create new one"""
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                return json.load(f)
        else:
            return {
                'project_name': self.project_name,
                'created': datetime.now().isoformat(),
                'version': '1.0.0',
                'files': {},
                'dependencies': {},
                'metadata': {}
            }
    
    def _save_registry(self):
        """Save registry to disk"""
        self.registry['last_modified'] = datetime.now().isoformat()
        with open(self.registry_file, 'w') as f:
            json.dump(self.registry, f, indent=2, default=str)
    
    def _generate_file_hash(self, data: Any) -> str:
        """Generate hash for data object"""
        if isinstance(data, pd.DataFrame):
            return hashlib.md5(pd.util.hash_pandas_object(data).values.tobytes()).hexdigest()[:8]
        elif isinstance(data, (dict, list)):
            return hashlib.md5(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()[:8]
        else:
            return hashlib.md5(str(data).encode()).hexdigest()[:8]
    
    def register_file(self, 
                      file_id: str,
                      data: Any,
                      file_type: str,
                      description: str = "",
                      dependencies: List[str] = None,
                      metadata: Dict = None) -> Path:
        """
        Register a file in the registry with metadata and dependencies
        
        Args:
            file_id: Unique identifier for the file
            data: Data to save (DataFrame, dict, etc.)
            file_type: Type of file (csv, bed, gff, json, etc.)
            description: Human-readable description
            dependencies: List of file_ids this file depends on
            metadata: Additional metadata dictionary
            
        Returns:
            Path to saved file
        """
        if dependencies is None:
            dependencies = []
        if metadata is None:
            metadata = {}
            
        # Generate file hash for integrity checking
        file_hash = self._generate_file_hash(data)
        
        # Determine file path based on type
        if file_type in ['bed', 'gff', 'json', 'fasta']:
            file_path = self.directories[file_type] / f"{file_id}.{file_type}"
        else:
            file_path = self.directories['data'] / f"{file_id}.{file_type}"
        
        # Save the actual file
        self._save_file(data, file_path, file_type)
        
        # Register in registry
        self.registry['files'][file_id] = {
            'path': str(file_path.relative_to(self.base_dir)),
            'type': file_type,
            'description': description,
            'hash': file_hash,
            'created': datetime.now().isoformat(),
            'size_bytes': file_path.stat().st_size if file_path.exists() else 0,
            'dependencies': dependencies,
            'metadata': metadata
        }
        
        # Update dependencies
        if dependencies:
            self.registry['dependencies'][file_id] = dependencies
        
        self._save_registry()
        logger.info(f"Registered file: {file_id} -> {file_path}")
        
        return file_path
    
    def _save_file(self, data: Any, file_path: Path, file_type: str):
        """Save data to appropriate file format"""
        if file_type == 'csv':
            if isinstance(data, pd.DataFrame):
                data.to_csv(file_path, index=False)
            else:
                raise ValueError(f"CSV format requires DataFrame, got {type(data)}")
        
        elif file_type == 'json':
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2, default=str)
        
        elif file_type == 'bed':
            self._save_bed_format(data, file_path)
        
        elif file_type == 'gff':
            self._save_gff_format(data, file_path)
        
        elif file_type == 'fasta':
            self._save_fasta_format(data, file_path)
        
        elif file_type == 'pkl':
            with open(file_path, 'wb') as f:
                pickle.dump(data, f)
        
        else:
            raise ValueError(f"Unsupported file type: {file_type}")
    
    def _save_bed_format(self, data: pd.DataFrame, file_path: Path):
        """Save DataFrame as BED format"""
        # Ensure required BED columns
        required_cols = ['chrom', 'start', 'end']
        if not all(col in data.columns for col in required_cols):
            raise ValueError(f"BED format requires columns: {required_cols}")
        
        # Standard BED format with optional additional columns
        bed_data = data.copy()
        
        # Ensure proper BED format (0-based coordinates)
        if 'start' in bed_data.columns:
            bed_data['start'] = bed_data['start'].astype(int)
        if 'end' in bed_data.columns:
            bed_data['end'] = bed_data['end'].astype(int)
        
        bed_data.to_csv(file_path, sep='\t', header=False, index=False)
    
    def _save_gff_format(self, data: pd.DataFrame, file_path: Path):
        """Save DataFrame as GFF3 format"""
        with open(file_path, 'w') as f:
            f.write("##gff-version 3\n")
            
            for _, row in data.iterrows():
                # Standard GFF3 fields
                fields = [
                    str(row.get('seqid', '.')),
                    str(row.get('source', 'genome_inversion_analyzer')),
                    str(row.get('type', 'feature')),
                    str(int(row.get('start', 1))),
                    str(int(row.get('end', 1))),
                    str(row.get('score', '.')),
                    str(row.get('strand', '.')),
                    str(row.get('phase', '.')),
                    str(row.get('attributes', '.'))
                ]
                f.write('\t'.join(fields) + '\n')
    
    def _save_fasta_format(self, data: Dict[str, str], file_path: Path):
        """Save sequences as FASTA format"""
        with open(file_path, 'w') as f:
            for seq_id, sequence in data.items():
                f.write(f">{seq_id}\n{sequence}\n")
    
    def get_file_path(self, file_id: str) -> Optional[Path]:
        """Get full path for registered file"""
        if file_id in self.registry['files']:
            rel_path = self.registry['files'][file_id]['path']
            return self.base_dir / rel_path
        return None
    
    def get_file_info(self, file_id: str) -> Optional[Dict]:
        """Get complete file information"""
        return self.registry['files'].get(file_id)
    
    def list_files(self, file_type: str = None) -> List[str]:
        """List all registered files, optionally filtered by type"""
        if file_type:
            return [fid for fid, info in self.registry['files'].items() 
                   if info['type'] == file_type]
        return list(self.registry['files'].keys())
    
    def get_dependencies(self, file_id: str) -> List[str]:
        """Get dependencies for a file"""
        return self.registry['dependencies'].get(file_id, [])
    
    def get_dependents(self, file_id: str) -> List[str]:
        """Get files that depend on this file"""
        dependents = []
        for fid, deps in self.registry['dependencies'].items():
            if file_id in deps:
                dependents.append(fid)
        return dependents
    
    def verify_integrity(self, file_id: str) -> bool:
        """Verify file integrity using stored hash"""
        file_info = self.get_file_info(file_id)
        if not file_info:
            return False
        
        file_path = self.get_file_path(file_id)
        if not file_path.exists():
            return False
        
        # Re-read and hash the file
        try:
            if file_info['type'] == 'csv':
                data = pd.read_csv(file_path)
            elif file_info['type'] == 'json':
                with open(file_path, 'r') as f:
                    data = json.load(f)
            elif file_info['type'] == 'pkl':
                with open(file_path, 'rb') as f:
                    data = pickle.load(f)
            else:
                # For text files, just check if they exist and are readable
                with open(file_path, 'r') as f:
                    f.read(1)  # Read one character to test
                return True
            
            current_hash = self._generate_file_hash(data)
            return current_hash == file_info['hash']
        
        except Exception as e:
            logger.warning(f"Failed to verify integrity for {file_id}: {e}")
            return False
    
    def export_manifest(self) -> Path:
        """Export complete file manifest"""
        manifest_path = self.base_dir / "file_manifest.json"
        
        manifest = {
            'project': self.registry['project_name'],
            'generated': datetime.now().isoformat(),
            'total_files': len(self.registry['files']),
            'files': []
        }
        
        for file_id, info in self.registry['files'].items():
            file_manifest = {
                'id': file_id,
                'path': info['path'],
                'type': info['type'],
                'description': info['description'],
                'size_mb': round(info['size_bytes'] / 1024 / 1024, 2),
                'dependencies': self.get_dependencies(file_id),
                'dependents': self.get_dependents(file_id)
            }
            manifest['files'].append(file_manifest)
        
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        logger.info(f"Exported manifest: {manifest_path}")
        return manifest_path
    
    def cleanup_orphaned_files(self, dry_run: bool = True) -> List[str]:
        """Remove files not referenced in registry"""
        orphaned = []
        
        for directory in self.directories.values():
            if directory.exists():
                for file_path in directory.iterdir():
                    if file_path.is_file() and file_path.name not in ['registry.json', 'file_manifest.json']:
                        # Check if file is registered
                        is_registered = any(
                            self.base_dir / info['path'] == file_path 
                            for info in self.registry['files'].values()
                        )
                        
                        if not is_registered:
                            orphaned.append(str(file_path))
                            if not dry_run:
                                file_path.unlink()
        
        if orphaned:
            action = "Would remove" if dry_run else "Removed"
            logger.info(f"{action} {len(orphaned)} orphaned files")
        
        return orphaned


class AnalysisResultsExporter:
    """
    Specialized exporter for genome inversion analysis results
    """
    
    def __init__(self, registry: FileRegistry):
        self.registry = registry
    
    def export_busco_coordinates(self, ortholog_df: pd.DataFrame) -> Path:
        """Export BUSCO coordinates as BED format"""
        bed_data = []
        
        for _, row in ortholog_df.iterrows():
            # First genome
            bed_data.append({
                'chrom': row['first_chr'],
                'start': row['first_start'] - 1,  # Convert to 0-based
                'end': row['first_end'],
                'name': f"{row['busco_id']}_first",
                'score': int(row.get('similarity', 0) * 1000),
                'strand': row['first_strand']
            })
            
            # Second genome
            bed_data.append({
                'chrom': row['second_chr'],
                'start': row['second_start'] - 1,
                'end': row['second_end'],
                'name': f"{row['busco_id']}_second",
                'score': int(row.get('similarity', 0) * 1000),
                'strand': row['second_strand']
            })
        
        bed_df = pd.DataFrame(bed_data)
        return self.registry.register_file(
            'busco_coordinates',
            bed_df,
            'bed',
            'BUSCO gene coordinates for genome browser visualization',
            dependencies=['ortholog_analysis'],
            metadata={'format': 'BED6', 'coordinate_system': '0-based'}
        )
    
    def export_synteny_blocks(self, synteny_df: pd.DataFrame) -> Path:
        """Export synteny blocks as BED format"""
        # This would need additional logic to convert synteny blocks to BED intervals
        # For now, export as JSON with proper structure
        synteny_data = {
            'synteny_blocks': synteny_df.to_dict('records'),
            'metadata': {
                'total_blocks': len(synteny_df),
                'format': 'synteny_json_v1'
            }
        }
        
        return self.registry.register_file(
            'synteny_blocks',
            synteny_data,
            'json',
            'Synteny block definitions with correlation and strand data',
            dependencies=['ortholog_analysis'],
            metadata={'format': 'synteny_json_v1'}
        )
    
    def export_inversion_regions(self, inversion_df: pd.DataFrame) -> Path:
        """Export inversion regions as BEDPE format for structural variants"""
        bedpe_data = []
        
        for _, row in inversion_df.iterrows():
            bedpe_data.append({
                'chrom1': row['first_chr'],
                'start1': 0,  # Would need actual breakpoint coordinates
                'end1': 1,
                'chrom2': row['second_chr'],
                'start2': 0,
                'end2': 1,
                'name': f"inversion_{row.get('start_gene', 'unknown')}",
                'score': int(row.get('confidence', 0) * 1000),
                'strand1': '+',
                'strand2': '-',
                'size_genes': row.get('size_genes', 0),
                'inversion_type': row.get('inversion_type', 'unknown')
            })
        
        bedpe_df = pd.DataFrame(bedpe_data)
        
        # Save as tab-separated BEDPE format manually (not through _save_file)
        bedpe_path = self.registry.directories['bed'] / 'inversions.bedpe'
        bedpe_df.to_csv(bedpe_path, sep='\t', index=False, header=True)
        
        # Register as JSON to avoid BED format validation issues
        return self.registry.register_file(
            'inversion_bedpe',
            bedpe_data,  # Store as list for JSON
            'json',
            'Inversion regions in BEDPE format for structural variant analysis',
            dependencies=['inversion_analysis'],
            metadata={'format': 'BEDPE', 'svtype': 'INV', 'bedpe_file': str(bedpe_path)}
        )
    
    def export_analysis_summary(self, results_dict: Dict) -> Path:
        """Export complete analysis summary as JSON"""
        summary = {
            'analysis_metadata': {
                'timestamp': datetime.now().isoformat(),
                'version': '1.0.0',
                'total_orthologs': len(results_dict.get('ortholog_df', [])),
                'total_synteny_blocks': len(results_dict.get('synteny_df', [])),
                'total_inversions': len(results_dict.get('inversion_df', [])),
                'total_rearrangements': len(results_dict.get('rearrangement_df', []))
            },
            'quality_metrics': {
                'first_genome': results_dict.get('first_quality', {}),
                'second_genome': results_dict.get('second_quality', {})
            }
        }
        
        return self.registry.register_file(
            'analysis_summary',
            summary,
            'json',
            'Complete analysis summary with metadata and quality metrics',
            dependencies=[],
            metadata={'format': 'analysis_summary_v1'}
        )