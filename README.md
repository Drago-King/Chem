# Chemistry Education Telegram Bot (@IUPAC_name_bot)

A comprehensive chemistry education Telegram bot featuring IUPAC name conversion, periodic table exploration, reaction animations, and an extensive compound database.

## Features

### 🧪 Chemical Conversion
- Convert IUPAC names to molecular structures
- Generate 2D and 3D molecular visualizations
- Calculate molecular properties (mass, bonds, LogP, etc.)
- Support for 119+ compounds across 20+ categories

### 🔬 Interactive Periodic Table
- High-quality 4K periodic table visualization
- Detailed element information with orbital diagrams
- Interactive element lookup with facts and properties
- Complete data for all 118 elements

### ⚗️ Chemical Reaction Animations
- 37 animated chemical reactions
- Universal GIF format for all devices
- Combustion, synthesis, decomposition, and biological processes
- Visual molecular transformation animations

### 📚 Compound Database
- 119 compounds across 20+ categories
- Enhanced molecular formula rendering (C₉H₈O₄ style)
- Pharmaceuticals, organic compounds, nucleotides, and more
- Searchable by name, formula, or category

## Commands

- `/start` - Welcome message and bot overview
- `/examples` - Show example compounds with quick buttons
- `/table` - Display interactive periodic table
- `/ask [element]` - Get detailed element information
- `/reactions` - View all available reaction animations
- `/animate [reaction]` - Create animated GIF of chemical reaction
- `/compounds [search]` - Search compound database
- Send any chemical name for instant conversion

## Installation

### Requirements
- Python 3.8+
- Telegram Bot Token
- Required Python packages (see requirements.txt)

### Setup
1. Clone this repository
2. Install dependencies: `pip install -r requirements.txt`
3. Set environment variables (see .env.example)
4. Run the bot: `python telegram_bot_simple.py`

### Environment Variables
Create a `.env` file with:
```
TELEGRAM_BOT_TOKEN=your_bot_token_here
```

## Project Structure

```
chemistry_telegram_bot/
├── telegram_bot_simple.py          # Main bot implementation
├── chemical_converter.py           # IUPAC name conversion
├── chemical_reaction_animator.py   # Animation generation
├── compound_visualizer.py          # Compound database visualization
├── common_compounds_database.py    # Compound data
├── periodic_table_visualizer.py   # Periodic table rendering
├── periodic_table_data.py          # Element data
├── molecular_3d_renderer.py        # 3D structure rendering
├── requirements.txt                # Python dependencies
├── .env.example                    # Environment variables template
└── assets/                         # Image assets
    └── periodic_table-4k.png       # High-quality periodic table
```

## Features Overview

### Chemical Name Conversion
- Supports IUPAC nomenclature
- Multiple data sources (PubChem, ChEBI, local database)
- Molecular structure visualization
- Property calculations

### Animation System
- 37 different chemical reactions
- Standard GIF format for universal compatibility
- Multiple fallback options for reliability
- Progress feedback during generation

### Compound Database
- 119+ compounds with full chemical data
- Categories: Pharmaceuticals, Organic Acids, Alcohols, etc.
- Enhanced molecular formula rendering
- Visual and text-based compound information

### Periodic Table
- Interactive element lookup
- Orbital configuration diagrams
- Element facts and properties
- High-quality 4K visualization

## Development

### Adding New Compounds
Edit `common_compounds_database.py` to add new compounds:
```python
'compound_id': {
    'iupac_name': 'IUPAC name',
    'common_name': 'Common name',
    'formula': 'C₆H₁₂O₆',
    'smiles': 'SMILES notation',
    'category': 'Category',
    'uses': 'Applications',
    'molar_mass': 180.16
}
```

### Adding New Reactions
Edit `chemical_reaction_animator.py` to add new reactions:
```python
'reaction_name': {
    'equation': 'Chemical equation',
    'description': 'Description',
    'energy_type': 'exothermic/endothermic'
}
```

## Deployment

### Local Development
1. Install dependencies
2. Set environment variables
3. Run: `python telegram_bot_simple.py`

### Production Deployment
- Supports continuous deployment
- Environment variable configuration
- Automatic restart on failure
- 24/7 operation capability

## License

This project is open source and available under the MIT License.

## Support

For questions, issues, or contributions, please open an issue on GitHub.

## Bot Commands Quick Reference

| Command | Description | Example |
|---------|-------------|---------|
| `/start` | Bot introduction | `/start` |
| `/examples` | Example compounds | `/examples` |
| `/table` | Periodic table | `/table` |
| `/ask` | Element info | `/ask hydrogen` |
| `/reactions` | Available animations | `/reactions` |
| `/animate` | Create animation | `/animate combustion_methane` |
| `/compounds` | Search compounds | `/compounds caffeine` |

## Technical Details

- **Language**: Python 3.8+
- **Framework**: Custom Telegram Bot API implementation
- **Chemistry Libraries**: RDKit, Matplotlib, Pillow
- **Image Processing**: PIL/Pillow for visualizations
- **Animation**: Matplotlib + PIL for GIF generation
- **Data Sources**: PubChem, ChEBI, local database

## Recent Updates

- Enhanced molecular formula rendering with clear subscript numbers
- Expanded compound database to 119+ compounds
- Fixed animation compatibility issues for universal device support
- Added comprehensive error handling and fallback options
- Improved visual quality for all compound displays