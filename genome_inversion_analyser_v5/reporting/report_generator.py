# =============================================================================
# Report Generator (reporting/report_generator.py)
# =============================================================================

"""
Comprehensive report generation system for genome analysis results.
Creates formatted text reports, HTML summaries, and statistical analyses.
"""

from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime
import json

from ..logger import get_logger

logger = get_logger()

class ReportGenerator:
    """
    Comprehensive report generator for genome analysis results.
    Creates multiple report formats including text, HTML, and JSON.
    """
    
    def __init__(self, config):
        """
        Initialize report generator.
        
        Args:
            config: Configuration object with reporting parameters
        """
        self.config = config
        self.output_dir = Path(config.get('output_directory', 'output'))
        self.create_html_report = config.get('create_html_report', True)
        self.create_json_report = config.get('create_json_report', True)
        self.include_detailed_stats = config.get('include_detailed_statistics', True)
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Report generator initialized")
        logger.info(f"  Output directory: {self.output_dir}")
        logger.info(f"  HTML report: {'enabled' if self.create_html_report else 'disabled'}")
        logger.info(f"  JSON report: {'enabled' if self.create_json_report else 'disabled'}")
    
    def generate_comprehensive_report(self, results: Dict[str, Any],
                                    plot_paths: Dict[str, str] = None) -> Dict[str, str]:
        """
        Generate comprehensive analysis report in multiple formats.
        
        Args:
            results: Complete analysis results dictionary
            plot_paths: Optional dictionary with paths to generated plots
            
        Returns:
            Dictionary with paths to generated report files
        """
        logger.info("Generating comprehensive analysis report...")
        
        report_paths = {}
        
        # 1. Generate text report
        text_report_path = self._generate_text_report(results, plot_paths)
        if text_report_path:
            report_paths['text_report'] = text_report_path
        
        # 2. Generate HTML report
        if self.create_html_report:
            html_report_path = self._generate_html_report(results, plot_paths)
            if html_report_path:
                report_paths['html_report'] = html_report_path
        
        # 3. Generate JSON report
        if self.create_json_report:
            json_report_path = self._generate_json_report(results)
            if json_report_path:
                report_paths['json_report'] = json_report_path
        
        # 4. Generate summary statistics file
        stats_path = self._generate_statistics_file(results)
        if stats_path:
            report_paths['statistics_file'] = stats_path
        
        logger.info(f"  Report generation completed: {len(report_paths)} files created")
        return report_paths
    
    def _generate_text_report(self, results: Dict[str, Any], 
                            plot_paths: Dict[str, str] = None) -> str:
        """Generate comprehensive text report."""
        report_path = str(self.output_dir / "genome_analysis_report.txt")
        
        try:
            with open(report_path, 'w') as f:
                # Header
                f.write("=" * 80 + "\n")
                f.write("GENOME COMPARISON ANALYSIS REPORT\n")
                f.write("=" * 80 + "\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                # Executive Summary
                self._write_executive_summary(f, results)
                
                # Assembly Quality Assessment
                self._write_quality_assessment(f, results)
                
                # Synteny Analysis
                self._write_synteny_analysis(f, results)
                
                # Inversion Analysis
                self._write_inversion_analysis(f, results)
                
                # Rearrangement Analysis
                self._write_rearrangement_analysis(f, results)
                
                # Statistical Validation
                self._write_statistical_validation(f, results)
                
                # Generated Files Summary
                if plot_paths:
                    self._write_files_summary(f, plot_paths)
                
                # Configuration Summary
                self._write_configuration_summary(f)
                
            logger.info(f"  Text report saved: {report_path}")
            return report_path
            
        except Exception as e:
            logger.error(f"Failed to generate text report: {e}")
            return None
    
    def _write_executive_summary(self, f, results: Dict[str, Any]):
        """Write executive summary section."""
        f.write("EXECUTIVE SUMMARY\n")
        f.write("-" * 40 + "\n\n")
        
        stats = results.get('statistics', {}).get('summary', {})
        
        # Main findings
        synteny_blocks = stats.get('total_synteny_blocks', 0)
        inversions = stats.get('total_inversions', 0)
        rearrangements = stats.get('total_rearrangements', 0)
        
        f.write(f"• Synteny blocks identified: {synteny_blocks}\n")
        f.write(f"• Inversions detected: {inversions}\n")
        f.write(f"• Chromosomal rearrangements found: {rearrangements}\n\n")
        
        # Quality assessment
        if 'first_quality' in results and 'second_quality' in results:
            q1 = results['first_quality']
            q2 = results['second_quality']
            f.write(f"• First genome quality: {q1['quality_class']} (score: {q1['quality_score']:.3f})\n")
            f.write(f"• Second genome quality: {q2['quality_class']} (score: {q2['quality_score']:.3f})\n\n")
        
        # Overall assessment
        if inversions > 0 or rearrangements > 0:
            f.write("The analysis reveals significant structural differences between the genomes,\n")
            f.write("including inversions and/or chromosomal rearrangements that suggest\n")
            f.write("evolutionary divergence or assembly differences.\n\n")
        else:
            f.write("The genomes show high structural similarity with minimal evidence\n")
            f.write("of major rearrangements or inversions.\n\n")
    
    def _write_quality_assessment(self, f, results: Dict[str, Any]):
        """Write assembly quality assessment section."""
        f.write("ASSEMBLY QUALITY ASSESSMENT\n")
        f.write("-" * 40 + "\n\n")
        
        if 'first_quality' in results:
            self._write_individual_quality(f, results['first_quality'], "First Genome")
        
        if 'second_quality' in results:
            self._write_individual_quality(f, results['second_quality'], "Second Genome")
        
        if 'quality_comparison' in results:
            comp = results['quality_comparison']
            f.write("Quality Comparison:\n")
            f.write(f"• Quality difference: {comp['quality_difference']:.3f}\n")
            for suggestion in comp.get('suggestions', []):
                f.write(f"• {suggestion}\n")
            f.write("\n")
    
    def _write_individual_quality(self, f, quality: Dict[str, Any], genome_name: str):
        """Write quality assessment for individual genome."""
        f.write(f"{genome_name} Quality Assessment:\n")
        f.write(f"• Overall quality: {quality['quality_class']} (score: {quality['quality_score']:.3f})\n")
        
        metrics = quality.get('metrics', {})
        if 'total_length' in metrics:
            f.write(f"• Total length: {metrics['total_length']:,} bp ({metrics['total_length']/1e6:.1f} Mbp)\n")
        if 'n_contigs' in metrics:
            f.write(f"• Number of contigs: {metrics['n_contigs']:,}\n")
        if 'n50' in metrics:
            f.write(f"• N50: {metrics['n50']:,} bp ({metrics['n50']/1e6:.1f} Mbp)\n")
        if 'busco_completeness' in metrics:
            f.write(f"• BUSCO completeness: {metrics['busco_completeness']:.1%}\n")
        if 'busco_duplication' in metrics:
            f.write(f"• BUSCO duplication: {metrics['busco_duplication']:.1%}\n")
        f.write("\n")
    
    def _write_synteny_analysis(self, f, results: Dict[str, Any]):
        """Write synteny analysis section."""
        f.write("SYNTENY ANALYSIS\n")
        f.write("-" * 40 + "\n\n")
        
        if 'synteny_df' in results and len(results['synteny_df']) > 0:
            synteny_df = results['synteny_df']
            f.write(f"Total synteny blocks: {len(synteny_df)}\n")
            
            if 'block_size' in synteny_df.columns:
                f.write(f"Average block size: {synteny_df['block_size'].mean():.1f} genes\n")
                f.write(f"Largest block: {synteny_df['block_size'].max()} genes\n")
                f.write(f"Smallest block: {synteny_df['block_size'].min()} genes\n")
            
            if 'synteny_type' in synteny_df.columns:
                type_counts = synteny_df['synteny_type'].value_counts()
                f.write("\nSynteny types:\n")
                for stype, count in type_counts.items():
                    f.write(f"• {stype}: {count} blocks ({count/len(synteny_df):.1%})\n")
            
            if 'confidence' in synteny_df.columns:
                f.write(f"\nAverage confidence: {synteny_df['confidence'].mean():.3f}\n")
                high_conf = (synteny_df['confidence'] >= 0.8).sum()
                f.write(f"High confidence blocks (≥0.8): {high_conf} ({high_conf/len(synteny_df):.1%})\n")
        else:
            f.write("No synteny blocks detected.\n")
        f.write("\n")
    
    def _write_inversion_analysis(self, f, results: Dict[str, Any]):
        """Write inversion analysis section."""
        f.write("INVERSION ANALYSIS\n")
        f.write("-" * 40 + "\n\n")
        
        if 'inversion_df' in results and len(results['inversion_df']) > 0:
            inversion_df = results['inversion_df']
            f.write(f"Total inversions detected: {len(inversion_df)}\n")
            
            if 'size_genes' in inversion_df.columns:
                f.write(f"Average inversion size: {inversion_df['size_genes'].mean():.1f} genes\n")
                f.write(f"Largest inversion: {inversion_df['size_genes'].max()} genes\n")
                total_inverted = inversion_df['size_genes'].sum()
                f.write(f"Total genes in inversions: {total_inverted}\n")
            
            if 'inversion_type' in inversion_df.columns:
                type_counts = inversion_df['inversion_type'].value_counts()
                f.write("\nInversion types:\n")
                for itype, count in type_counts.items():
                    f.write(f"• {itype.replace('_', ' ')}: {count} ({count/len(inversion_df):.1%})\n")
            
            if 'detection_method' in inversion_df.columns:
                method_counts = inversion_df['detection_method'].value_counts()
                f.write("\nDetection methods:\n")
                for method, count in method_counts.items():
                    f.write(f"• {method.replace('_', ' ')}: {count} inversions\n")
            
            if 'confidence' in inversion_df.columns:
                f.write(f"\nAverage confidence: {inversion_df['confidence'].mean():.3f}\n")
                high_conf = (inversion_df['confidence'] >= 0.8).sum()
                f.write(f"High confidence inversions (≥0.8): {high_conf} ({high_conf/len(inversion_df):.1%})\n")
        else:
            f.write("No inversions detected.\n")
        f.write("\n")
    
    def _write_rearrangement_analysis(self, f, results: Dict[str, Any]):
        """Write rearrangement analysis section."""
        f.write("CHROMOSOMAL REARRANGEMENT ANALYSIS\n")
        f.write("-" * 40 + "\n\n")
        
        if 'rearrangement_df' in results and len(results['rearrangement_df']) > 0:
            rearr_df = results['rearrangement_df']
            f.write(f"Total rearrangements detected: {len(rearr_df)}\n")
            
            if 'type' in rearr_df.columns:
                type_counts = rearr_df['type'].value_counts()
                f.write("\nRearrangement types:\n")
                for rtype, count in type_counts.items():
                    f.write(f"• {rtype.replace('_', ' ').title()}: {count}\n")
                    
                    # Add specific details for each type
                    type_data = rearr_df[rearr_df['type'] == rtype]
                    if rtype == 'chromosome_split' and 'split_ratio' in type_data.columns:
                        avg_ratio = type_data['split_ratio'].mean()
                        f.write(f"  - Average split ratio: {avg_ratio:.1f}\n")
                    elif rtype == 'chromosome_fusion' and 'fusion_ratio' in type_data.columns:
                        avg_ratio = type_data['fusion_ratio'].mean()
                        f.write(f"  - Average fusion ratio: {avg_ratio:.1f}\n")
                    elif rtype == 'reciprocal_translocation' and 'reciprocity_balance' in type_data.columns:
                        avg_balance = type_data['reciprocity_balance'].mean()
                        f.write(f"  - Average reciprocity balance: {avg_balance:.3f}\n")
            
            if 'confidence' in rearr_df.columns:
                f.write(f"\nAverage confidence: {rearr_df['confidence'].mean():.3f}\n")
                high_conf = (rearr_df['confidence'] >= 0.8).sum()
                f.write(f"High confidence rearrangements (≥0.8): {high_conf} ({high_conf/len(rearr_df):.1%})\n")
        else:
            f.write("No chromosomal rearrangements detected.\n")
        f.write("\n")
    
    def _write_statistical_validation(self, f, results: Dict[str, Any]):
        """Write statistical validation section."""
        f.write("STATISTICAL VALIDATION\n")
        f.write("-" * 40 + "\n\n")
        
        # Synteny validation
        if 'synteny_validation' in results:
            val = results['synteny_validation']
            overall = val.get('overall_validation', {})
            f.write(f"Synteny validation: {'PASSED' if overall.get('validated', False) else 'FAILED'}\n")
            f.write(f"Confidence level: {overall.get('confidence_level', 'N/A')}\n")
            
            if 'correlation_validation' in val:
                corr_val = val['correlation_validation']
                if 'mean_correlation' in corr_val:
                    f.write(f"Mean position correlation: {corr_val['mean_correlation']:.3f}\n")
                if 'p_value' in corr_val:
                    f.write(f"Correlation significance (p-value): {corr_val['p_value']:.3e}\n")
        
        # Inversion validation
        if 'inversion_validation' in results:
            val = results['inversion_validation']
            overall = val.get('overall_validation', {})
            f.write(f"Inversion validation: {'PASSED' if overall.get('validated', False) else 'FAILED'}\n")
        
        if 'synteny_validation' not in results and 'inversion_validation' not in results:
            f.write("Statistical validation was not performed.\n")
        
        f.write("\n")
    
    def _write_files_summary(self, f, plot_paths: Dict[str, str]):
        """Write generated files summary."""
        f.write("GENERATED FILES\n")
        f.write("-" * 40 + "\n\n")
        
        f.write("Visualization files:\n")
        for plot_type, plot_path in plot_paths.items():
            if not plot_type.startswith('chr_plot_'):
                plot_name = plot_type.replace('_', ' ').title()
                f.write(f"• {plot_name}: {Path(plot_path).name}\n")
        
        # Count individual chromosome plots
        chr_plots = [k for k in plot_paths.keys() if k.startswith('chr_plot_')]
        if chr_plots:
            f.write(f"• Individual chromosome plots: {len(chr_plots)} files\n")
        
        f.write("\n")
    
    def _write_configuration_summary(self, f):
        """Write analysis configuration summary."""
        f.write("ANALYSIS CONFIGURATION\n")
        f.write("-" * 40 + "\n\n")
        
        # Key parameters
        config_items = [
            ('Minimum synteny block size', 'base_min_synteny_block_size', 3),
            ('Synteny correlation threshold', 'base_synteny_correlation_threshold', 0.8),
            ('Minimum inversion size', 'base_min_inversion_size', 2),
            ('BUSCO similarity threshold', 'busco_similarity_threshold', 0.8),
            ('Statistical validation', 'enable_statistical_validation', False),
            ('Single gene inversions', 'enable_single_gene_inversions', True),
            ('Paralog detection', 'enable_paralog_detection', True)
        ]
        
        for name, key, default in config_items:
            value = self.config.get(key, default)
            f.write(f"• {name}: {value}\n")
        
        f.write(f"\nAnalysis performed with {len(self.config)} configuration parameters.\n")
        f.write("For complete configuration details, see the JSON report.\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")
    
    def _generate_html_report(self, results: Dict[str, Any], 
                            plot_paths: Dict[str, str] = None) -> str:
        """Generate HTML report with embedded plots."""
        report_path = str(self.output_dir / "genome_analysis_report.html")
        
        try:
            with open(report_path, 'w') as f:
                # HTML header
                f.write(self._get_html_header())
                
                # Title and summary
                f.write(self._get_html_summary(results))
                
                # Quality assessment
                f.write(self._get_html_quality_section(results))
                
                # Analysis results
                f.write(self._get_html_analysis_section(results))
                
                # Visualizations
                if plot_paths:
                    f.write(self._get_html_visualizations_section(plot_paths))
                
                # Statistical details
                if self.include_detailed_stats:
                    f.write(self._get_html_statistics_section(results))
                
                # HTML footer
                f.write(self._get_html_footer())
            
            logger.info(f"  HTML report saved: {report_path}")
            return report_path
            
        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
            return None
    
    def _get_html_header(self) -> str:
        """Get HTML header with CSS styling."""
        return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome Comparison Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; color: #333; }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }
        .header h1 { margin: 0; font-size: 2.5em; text-align: center; }
        .header .subtitle { text-align: center; font-size: 1.2em; margin-top: 10px; opacity: 0.9; }
        .section { background: #f8f9fa; padding: 25px; margin: 20px 0; border-radius: 8px; border-left: 5px solid #007bff; }
        .section h2 { color: #007bff; margin-top: 0; border-bottom: 2px solid #e9ecef; padding-bottom: 10px; }
        .section h3 { color: #495057; margin-top: 20px; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0; }
        .stat-card { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); border-left: 4px solid #28a745; }
        .stat-card h4 { margin: 0 0 10px 0; color: #28a745; font-size: 1.1em; }
        .stat-card .value { font-size: 2em; font-weight: bold; color: #333; }
        .stat-card .description { color: #666; font-size: 0.9em; margin-top: 5px; }
        .quality-indicator { display: inline-block; padding: 5px 15px; border-radius: 20px; font-weight: bold; color: white; }
        .quality-high { background-color: #28a745; }
        .quality-medium { background-color: #ffc107; color: #212529; }
        .quality-low { background-color: #dc3545; }
        .quality-fragmented { background-color: #6c757d; }
        .plot-container { text-align: center; margin: 20px 0; }
        .plot-container img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
        .validation-passed { color: #28a745; font-weight: bold; }
        .validation-failed { color: #dc3545; font-weight: bold; }
        table { width: 100%; border-collapse: collapse; margin: 15px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f8f9fa; font-weight: bold; color: #495057; }
        .footer { text-align: center; margin-top: 50px; padding: 20px; background-color: #f8f9fa; border-radius: 8px; color: #666; }
    </style>
</head>
<body>
"""
    
    def _get_html_summary(self, results: Dict[str, Any]) -> str:
        """Get HTML summary section."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        stats = results.get('statistics', {}).get('summary', {})
        
        html = f"""
    <div class="header">
        <h1>Genome Comparison Analysis Report</h1>
        <div class="subtitle">Generated on {timestamp}</div>
    </div>
    
    <div class="section">
        <h2>Executive Summary</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h4>Synteny Blocks</h4>
                <div class="value">{stats.get('total_synteny_blocks', 0)}</div>
                <div class="description">Conserved gene order regions</div>
            </div>
            <div class="stat-card">
                <h4>Inversions</h4>
                <div class="value">{stats.get('total_inversions', 0)}</div>
                <div class="description">Detected inversion events</div>
            </div>
            <div class="stat-card">
                <h4>Rearrangements</h4>
                <div class="value">{stats.get('total_rearrangements', 0)}</div>
                <div class="description">Chromosomal rearrangements</div>
            </div>
        """
        
        if 'average_quality_score' in stats:
            html += f"""
            <div class="stat-card">
                <h4>Average Quality</h4>
                <div class="value">{stats['average_quality_score']:.3f}</div>
                <div class="description">Assembly quality score</div>
            </div>
        """
        
        html += """
        </div>
    </div>
"""
        return html
    
    def _get_html_quality_section(self, results: Dict[str, Any]) -> str:
        """Get HTML quality assessment section."""
        html = """
    <div class="section">
        <h2>Assembly Quality Assessment</h2>
"""
        
        if 'first_quality' in results and 'second_quality' in results:
            q1 = results['first_quality']
            q2 = results['second_quality']
            
            html += f"""
        <div class="stats-grid">
            <div class="stat-card">
                <h4>First Genome</h4>
                <div class="value">{q1['quality_score']:.3f}</div>
                <div class="description">
                    <span class="quality-indicator quality-{q1['quality_class']}">{q1['quality_class'].title()}</span>
                </div>
            </div>
            <div class="stat-card">
                <h4>Second Genome</h4>
                <div class="value">{q2['quality_score']:.3f}</div>
                <div class="description">
                    <span class="quality-indicator quality-{q2['quality_class']}">{q2['quality_class'].title()}</span>
                </div>
            </div>
        </div>
"""
            
            # Detailed metrics table
            html += """
        <h3>Detailed Quality Metrics</h3>
        <table>
            <tr><th>Metric</th><th>First Genome</th><th>Second Genome</th></tr>
"""
            
            metrics1 = q1.get('metrics', {})
            metrics2 = q2.get('metrics', {})
            
            metric_pairs = [
                ('Total Length (Mbp)', 'total_length', lambda x: f"{x/1e6:.1f}" if x else "N/A"),
                ('Number of Contigs', 'n_contigs', lambda x: f"{x:,}" if x else "N/A"),
                ('N50 (Mbp)', 'n50', lambda x: f"{x/1e6:.1f}" if x else "N/A"),
                ('BUSCO Completeness (%)', 'busco_completeness', lambda x: f"{x*100:.1f}" if x else "N/A"),
                ('BUSCO Duplication (%)', 'busco_duplication', lambda x: f"{x*100:.1f}" if x else "N/A")
            ]
            
            for name, key, formatter in metric_pairs:
                val1 = formatter(metrics1.get(key, 0))
                val2 = formatter(metrics2.get(key, 0))
                html += f"            <tr><td>{name}</td><td>{val1}</td><td>{val2}</td></tr>\n"
            
            html += "        </table>\n"
        
        html += "    </div>\n"
        return html
    
    def _get_html_analysis_section(self, results: Dict[str, Any]) -> str:
        """Get HTML analysis results section."""
        html = """
    <div class="section">
        <h2>Analysis Results</h2>
"""
        
        # Synteny results
        if 'synteny_df' in results and len(results['synteny_df']) > 0:
            synteny_df = results['synteny_df']
            html += f"""
        <h3>Synteny Analysis</h3>
        <p><strong>{len(synteny_df)} synteny blocks</strong> identified with an average size of 
        <strong>{synteny_df['block_size'].mean():.1f} genes</strong>.</p>
"""
            
            if 'synteny_type' in synteny_df.columns:
                type_counts = synteny_df['synteny_type'].value_counts()
                html += "        <ul>\n"
                for stype, count in type_counts.items():
                    html += f"            <li>{stype.replace('_', ' ').title()}: {count} blocks ({count/len(synteny_df):.1%})</li>\n"
                html += "        </ul>\n"
        
        # Inversion results
        if 'inversion_df' in results and len(results['inversion_df']) > 0:
            inversion_df = results['inversion_df']
            html += f"""
        <h3>Inversion Analysis</h3>
        <p><strong>{len(inversion_df)} inversions</strong> detected with an average size of 
        <strong>{inversion_df['size_genes'].mean():.1f} genes</strong>.</p>
"""
            
            if 'inversion_type' in inversion_df.columns:
                type_counts = inversion