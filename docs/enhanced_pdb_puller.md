# Enhanced PDB Puller Documentation

## Overview

The Enhanced PDB Puller is a comprehensive data fetching system that retrieves extensive information about protein structures from the RCSB Protein Data Bank (PDB). Unlike basic PDB fetching that only downloads structure files, the enhanced puller collects:

- **Structure files** (PDB format, optionally mmCIF)
- **Comprehensive metadata** via RCSB PDB REST API
- **Experimental details** and validation metrics
- **Sequence information** (FASTA format)
- **Ligand and heteroatom information**
- **Publication data** and cross-references
- **Validation reports** and quality metrics

## Features

### 1. Comprehensive Data Collection

The enhanced puller fetches data from multiple RCSB PDB API endpoints:

- **Entry Information**: Basic structure details, title, authors, dates
- **Experimental Data**: Method, resolution, R-factors, space group, unit cell
- **Polymer Entities**: Protein chain information and descriptions
- **Non-polymer Entities**: Ligands, cofactors, and other molecules
- **Publications**: Associated research papers and citations
- **Validation Data**: Quality metrics and validation reports
- **Assembly Information**: Biological assemblies and quaternary structure

### 2. Organized Data Storage

Downloaded data is organized in a structured directory hierarchy:

```
data/pulled_structures/
├── structures/          # PDB and mmCIF structure files
├── metadata/           # JSON metadata files
├── sequences/          # FASTA sequence files
└── validation/         # Validation reports and quality data
```

### 3. Enhanced User Interface

#### Fetch Summary Dialog
When fetching a PDB structure, users see a comprehensive summary including:
- Downloaded file types
- Metadata sections retrieved
- Basic structure information (title, method, resolution)
- Number of protein chains and ligands

#### PDB Information Dialog (Tools → Show PDB Information)
Accessible via `Ctrl+I`, this dialog provides:
- **Summary Tab**: Key structure information and overview
- **Detailed Metadata Tab**: Complete API response data
- **Sequences Tab**: Full protein sequences for all chains
- **Downloaded Files Tab**: Information about local files

#### Enhanced Structural Analysis
The structural analysis module now displays:
- **Enhanced Structure Information**: Tabbed interface with comprehensive metadata
- **Experimental Details**: Method, resolution, R-factors, space group, unit cell
- **Entities & Ligands**: Detailed tables of protein chains and bound molecules
- **Publications**: Associated research papers and citations

## Usage

### Basic Usage

1. **Fetch PDB Structure**: Enter a PDB ID and click "Fetch"
   - Downloads structure file and comprehensive metadata
   - Shows summary dialog with fetched information
   - Stores data for future access

2. **View PDB Information**: Use `Tools → Show PDB Information` or `Ctrl+I`
   - Displays comprehensive information about current structure
   - Shows cached data if available

3. **Enhanced Analysis**: Use the Structural Analysis tab
   - Automatically uses enhanced metadata when available
   - Displays comprehensive experimental and publication information

### Advanced Features

#### Programmatic Access

```python
from src.core.enhanced_pdb_puller import EnhancedPDBPuller

# Initialize puller
puller = EnhancedPDBPuller("data/pulled_structures")

# Fetch comprehensive data
data = puller.fetch_comprehensive_pdb_data(
    "1CRN",
    include_validation=True,
    include_sequences=True,
    include_mmcif=False
)

# Access specific information
metadata = data['metadata']
sequences = data['sequences']
files = data['files']
```

#### Configuration Options

- **include_validation**: Download validation reports and quality metrics
- **include_sequences**: Download FASTA sequence files
- **include_mmcif**: Download mmCIF format files (more comprehensive than PDB)

## Data Sources

### RCSB PDB REST API Endpoints

The enhanced puller uses the following API endpoints:

- `https://data.rcsb.org/rest/v1/core/entry/{pdb_id}` - Basic entry information
- `https://data.rcsb.org/rest/v1/core/experimental_method/{pdb_id}` - Experimental details
- `https://data.rcsb.org/rest/v1/core/summary/{pdb_id}` - Structure summary
- `https://data.rcsb.org/rest/v1/core/pubmed/{pdb_id}` - Publication information
- `https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}` - Protein chain data
- `https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}` - Ligand information
- `https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}` - Assembly information
- `https://data.rcsb.org/rest/v1/core/validation/{pdb_id}` - Validation data

### File Downloads

- **PDB Files**: `https://files.rcsb.org/download/{pdb_id}.pdb`
- **mmCIF Files**: `https://files.rcsb.org/download/{pdb_id}.cif`
- **FASTA Sequences**: `https://files.rcsb.org/download/{pdb_id}.fasta`
- **Validation Reports**: `https://files.rcsb.org/download/{pdb_id}_validation.pdf`

## Error Handling

The enhanced puller includes robust error handling:

- **Network Timeouts**: Configurable timeout for API requests
- **Missing Data**: Graceful handling of unavailable endpoints
- **File Download Errors**: Continues operation if optional files fail
- **API Rate Limiting**: Respects RCSB PDB API guidelines

## Benefits

### For Researchers

1. **Complete Information**: Access to all available PDB metadata in one operation
2. **Publication Tracking**: Automatic retrieval of associated research papers
3. **Quality Assessment**: Built-in access to validation metrics
4. **Ligand Analysis**: Comprehensive information about bound molecules

### For Educators

1. **Rich Context**: Students see complete experimental context
2. **Historical Information**: Access to deposition and publication dates
3. **Method Comparison**: Easy comparison of experimental techniques
4. **Quality Discussion**: Built-in quality metrics for teaching

### For Developers

1. **Structured Data**: Well-organized JSON metadata for programmatic access
2. **Caching**: Automatic caching prevents redundant downloads
3. **Extensible**: Easy to add new data sources and endpoints
4. **Error Resilient**: Robust error handling for production use

## Technical Details

### Dependencies

- **requests**: HTTP client for API calls
- **biopython**: PDB file handling and parsing
- **PyQt5**: User interface components

### Performance

- **Parallel Downloads**: Multiple files downloaded concurrently
- **Caching**: Previously downloaded data is reused
- **Incremental Updates**: Only missing data is fetched
- **Memory Efficient**: Large files streamed to disk

### Data Format

Metadata is stored as JSON files with the following structure:

```json
{
  "pdb_id": "1CRN",
  "fetch_timestamp": "2025-01-20 10:30:00",
  "files": {
    "pdb": "/path/to/1CRN.pdb",
    "fasta": "/path/to/1CRN.fasta"
  },
  "metadata": {
    "entry": { ... },
    "experimental": { ... },
    "publications": [ ... ],
    "polymer_entities": [ ... ],
    "nonpolymer_entities": [ ... ]
  },
  "sequences": {
    "chains": [ ... ]
  },
  "validation": {
    "report": { ... }
  }
}
```

## Future Enhancements

Planned improvements include:

1. **AlphaFold Integration**: Fetch AlphaFold predictions for structures
2. **ChEMBL Integration**: Retrieve bioactivity data for ligands
3. **UniProt Mapping**: Automatic protein sequence annotation
4. **Structure Comparison**: Compare multiple structures automatically
5. **Batch Processing**: Download multiple structures simultaneously

## Troubleshooting

### Common Issues

1. **Network Connectivity**: Ensure internet access for API calls
2. **Disk Space**: Check available space for downloaded files
3. **API Availability**: Some older structures may have limited metadata
4. **Rate Limiting**: Avoid rapid successive requests to respect API limits

### Debug Information

Enable debug output by setting environment variable:
```bash
export PICOMOL_DEBUG=1
```

This provides detailed logging of API requests and responses.