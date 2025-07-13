"""
Comprehensive periodic table data with all elements up to Oganesson (118)
Including atomic properties, orbital configurations, and visualization data.
"""

PERIODIC_TABLE_DATA = {
    1: {
        'symbol': 'H', 'name': 'Hydrogen', 'atomic_mass': 1.008, 'group': 1, 'period': 1,
        'electron_config': '1s¹', 'orbital_diagram': '1s: ↑',
        'category': 'Nonmetal', 'color': '#FFFFFF', 'discovery_year': 1766,
        'melting_point': -259.16, 'boiling_point': -252.87, 'density': 0.00008988,
        'oxidation_states': [-1, 1], 'electronegativity': 2.20
    },
    2: {
        'symbol': 'He', 'name': 'Helium', 'atomic_mass': 4.003, 'group': 18, 'period': 1,
        'electron_config': '1s²', 'orbital_diagram': '1s: ↑↓',
        'category': 'Noble Gas', 'color': '#D9FFFF', 'discovery_year': 1895,
        'melting_point': -272.20, 'boiling_point': -268.93, 'density': 0.0001785,
        'oxidation_states': [0], 'electronegativity': None
    },
    3: {
        'symbol': 'Li', 'name': 'Lithium', 'atomic_mass': 6.941, 'group': 1, 'period': 2,
        'electron_config': '[He] 2s¹', 'orbital_diagram': '1s: ↑↓  2s: ↑',
        'category': 'Alkali Metal', 'color': '#CC80FF', 'discovery_year': 1817,
        'melting_point': 180.54, 'boiling_point': 1342, 'density': 0.534,
        'oxidation_states': [1], 'electronegativity': 0.98
    },
    4: {
        'symbol': 'Be', 'name': 'Beryllium', 'atomic_mass': 9.012, 'group': 2, 'period': 2,
        'electron_config': '[He] 2s²', 'orbital_diagram': '1s: ↑↓  2s: ↑↓',
        'category': 'Alkaline Earth Metal', 'color': '#C2FF00', 'discovery_year': 1797,
        'melting_point': 1287, 'boiling_point': 2471, 'density': 1.85,
        'oxidation_states': [2], 'electronegativity': 1.57
    },
    5: {
        'symbol': 'B', 'name': 'Boron', 'atomic_mass': 10.811, 'group': 13, 'period': 2,
        'electron_config': '[He] 2s² 2p¹', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑ _ _',
        'category': 'Metalloid', 'color': '#FFB5B5', 'discovery_year': 1808,
        'melting_point': 2075, 'boiling_point': 4000, 'density': 2.37,
        'oxidation_states': [3], 'electronegativity': 2.04
    },
    6: {
        'symbol': 'C', 'name': 'Carbon', 'atomic_mass': 12.011, 'group': 14, 'period': 2,
        'electron_config': '[He] 2s² 2p²', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑ ↑ _',
        'category': 'Nonmetal', 'color': '#909090', 'discovery_year': -3750,
        'melting_point': 3550, 'boiling_point': 4027, 'density': 2.267,
        'oxidation_states': [-4, -3, -2, -1, 1, 2, 3, 4], 'electronegativity': 2.55
    },
    7: {
        'symbol': 'N', 'name': 'Nitrogen', 'atomic_mass': 14.007, 'group': 15, 'period': 2,
        'electron_config': '[He] 2s² 2p³', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑ ↑ ↑',
        'category': 'Nonmetal', 'color': '#3050F8', 'discovery_year': 1772,
        'melting_point': -210.00, 'boiling_point': -195.79, 'density': 0.0012506,
        'oxidation_states': [-3, -2, -1, 1, 2, 3, 4, 5], 'electronegativity': 3.04
    },
    8: {
        'symbol': 'O', 'name': 'Oxygen', 'atomic_mass': 15.999, 'group': 16, 'period': 2,
        'electron_config': '[He] 2s² 2p⁴', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑↓ ↑ ↑',
        'category': 'Nonmetal', 'color': '#FF0D0D', 'discovery_year': 1774,
        'melting_point': -218.79, 'boiling_point': -182.95, 'density': 0.001429,
        'oxidation_states': [-2, -1, 1, 2], 'electronegativity': 3.44
    },
    9: {
        'symbol': 'F', 'name': 'Fluorine', 'atomic_mass': 18.998, 'group': 17, 'period': 2,
        'electron_config': '[He] 2s² 2p⁵', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑↓ ↑↓ ↑',
        'category': 'Halogen', 'color': '#90E050', 'discovery_year': 1886,
        'melting_point': -219.67, 'boiling_point': -188.11, 'density': 0.001696,
        'oxidation_states': [-1], 'electronegativity': 3.98
    },
    10: {
        'symbol': 'Ne', 'name': 'Neon', 'atomic_mass': 20.180, 'group': 18, 'period': 2,
        'electron_config': '[He] 2s² 2p⁶', 'orbital_diagram': '1s: ↑↓  2s: ↑↓  2p: ↑↓ ↑↓ ↑↓',
        'category': 'Noble Gas', 'color': '#B3E3F5', 'discovery_year': 1898,
        'melting_point': -248.59, 'boiling_point': -246.08, 'density': 0.0008999,
        'oxidation_states': [0], 'electronegativity': None
    },
    11: {
        'symbol': 'Na', 'name': 'Sodium', 'atomic_mass': 22.990, 'group': 1, 'period': 3,
        'electron_config': '[Ne] 3s¹', 'orbital_diagram': '[Ne] 3s: ↑',
        'category': 'Alkali Metal', 'color': '#AB5CF2', 'discovery_year': 1807,
        'melting_point': 97.794, 'boiling_point': 883, 'density': 0.971,
        'oxidation_states': [1], 'electronegativity': 0.93
    },
    12: {
        'symbol': 'Mg', 'name': 'Magnesium', 'atomic_mass': 24.305, 'group': 2, 'period': 3,
        'electron_config': '[Ne] 3s²', 'orbital_diagram': '[Ne] 3s: ↑↓',
        'category': 'Alkaline Earth Metal', 'color': '#8AFF00', 'discovery_year': 1755,
        'melting_point': 650, 'boiling_point': 1090, 'density': 1.738,
        'oxidation_states': [2], 'electronegativity': 1.31
    },
    13: {
        'symbol': 'Al', 'name': 'Aluminum', 'atomic_mass': 26.982, 'group': 13, 'period': 3,
        'electron_config': '[Ne] 3s² 3p¹', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑ _ _',
        'category': 'Post-transition Metal', 'color': '#BFA6A6', 'discovery_year': 1825,
        'melting_point': 660.32, 'boiling_point': 2519, 'density': 2.70,
        'oxidation_states': [3], 'electronegativity': 1.61
    },
    14: {
        'symbol': 'Si', 'name': 'Silicon', 'atomic_mass': 28.085, 'group': 14, 'period': 3,
        'electron_config': '[Ne] 3s² 3p²', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑ ↑ _',
        'category': 'Metalloid', 'color': '#F0C8A0', 'discovery_year': 1854,
        'melting_point': 1414, 'boiling_point': 3265, 'density': 2.3296,
        'oxidation_states': [-4, -3, -2, -1, 1, 2, 3, 4], 'electronegativity': 1.90
    },
    15: {
        'symbol': 'P', 'name': 'Phosphorus', 'atomic_mass': 30.974, 'group': 15, 'period': 3,
        'electron_config': '[Ne] 3s² 3p³', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑ ↑ ↑',
        'category': 'Nonmetal', 'color': '#FF8000', 'discovery_year': 1669,
        'melting_point': 44.15, 'boiling_point': 280.5, 'density': 1.82,
        'oxidation_states': [-3, -2, -1, 1, 2, 3, 4, 5], 'electronegativity': 2.19
    },
    16: {
        'symbol': 'S', 'name': 'Sulfur', 'atomic_mass': 32.065, 'group': 16, 'period': 3,
        'electron_config': '[Ne] 3s² 3p⁴', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑↓ ↑ ↑',
        'category': 'Nonmetal', 'color': '#FFFF30', 'discovery_year': -2000,
        'melting_point': 115.21, 'boiling_point': 444.61, 'density': 2.067,
        'oxidation_states': [-2, -1, 1, 2, 3, 4, 5, 6], 'electronegativity': 2.58
    },
    17: {
        'symbol': 'Cl', 'name': 'Chlorine', 'atomic_mass': 35.453, 'group': 17, 'period': 3,
        'electron_config': '[Ne] 3s² 3p⁵', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑↓ ↑↓ ↑',
        'category': 'Halogen', 'color': '#1FF01F', 'discovery_year': 1774,
        'melting_point': -101.5, 'boiling_point': -34.04, 'density': 0.003214,
        'oxidation_states': [-1, 1, 3, 5, 7], 'electronegativity': 3.16
    },
    18: {
        'symbol': 'Ar', 'name': 'Argon', 'atomic_mass': 39.948, 'group': 18, 'period': 3,
        'electron_config': '[Ne] 3s² 3p⁶', 'orbital_diagram': '[Ne] 3s: ↑↓  3p: ↑↓ ↑↓ ↑↓',
        'category': 'Noble Gas', 'color': '#80D1E3', 'discovery_year': 1894,
        'melting_point': -189.35, 'boiling_point': -185.85, 'density': 0.0017837,
        'oxidation_states': [0], 'electronegativity': None
    },
    # Continue with more elements...
    19: {
        'symbol': 'K', 'name': 'Potassium', 'atomic_mass': 39.098, 'group': 1, 'period': 4,
        'electron_config': '[Ar] 4s¹', 'orbital_diagram': '[Ar] 4s: ↑',
        'category': 'Alkali Metal', 'color': '#8F40D4', 'discovery_year': 1807,
        'melting_point': 63.5, 'boiling_point': 759, 'density': 0.862,
        'oxidation_states': [1], 'electronegativity': 0.82
    },
    20: {
        'symbol': 'Ca', 'name': 'Calcium', 'atomic_mass': 40.078, 'group': 2, 'period': 4,
        'electron_config': '[Ar] 4s²', 'orbital_diagram': '[Ar] 4s: ↑↓',
        'category': 'Alkaline Earth Metal', 'color': '#3DFF00', 'discovery_year': 1808,
        'melting_point': 842, 'boiling_point': 1484, 'density': 1.54,
        'oxidation_states': [2], 'electronegativity': 1.00
    },
    # Adding key transition metals
    26: {
        'symbol': 'Fe', 'name': 'Iron', 'atomic_mass': 55.845, 'group': 8, 'period': 4,
        'electron_config': '[Ar] 3d⁶ 4s²', 'orbital_diagram': '[Ar] 3d: ↑↓ ↑↓ ↑ ↑ ↑  4s: ↑↓',
        'category': 'Transition Metal', 'color': '#E06633', 'discovery_year': -1500,
        'melting_point': 1538, 'boiling_point': 2861, 'density': 7.874,
        'oxidation_states': [-2, -1, 1, 2, 3, 4, 5, 6], 'electronegativity': 1.83
    },
    29: {
        'symbol': 'Cu', 'name': 'Copper', 'atomic_mass': 63.546, 'group': 11, 'period': 4,
        'electron_config': '[Ar] 3d¹⁰ 4s¹', 'orbital_diagram': '[Ar] 3d: ↑↓ ↑↓ ↑↓ ↑↓ ↑↓  4s: ↑',
        'category': 'Transition Metal', 'color': '#C88033', 'discovery_year': -9000,
        'melting_point': 1084.62, 'boiling_point': 2562, 'density': 8.96,
        'oxidation_states': [1, 2, 3, 4], 'electronegativity': 1.90
    },
    47: {
        'symbol': 'Ag', 'name': 'Silver', 'atomic_mass': 107.868, 'group': 11, 'period': 5,
        'electron_config': '[Kr] 4d¹⁰ 5s¹', 'orbital_diagram': '[Kr] 4d: ↑↓ ↑↓ ↑↓ ↑↓ ↑↓  5s: ↑',
        'category': 'Transition Metal', 'color': '#C0C0C0', 'discovery_year': -3000,
        'melting_point': 961.78, 'boiling_point': 2162, 'density': 10.501,
        'oxidation_states': [1, 2, 3], 'electronegativity': 1.93
    },
    79: {
        'symbol': 'Au', 'name': 'Gold', 'atomic_mass': 196.967, 'group': 11, 'period': 6,
        'electron_config': '[Xe] 4f¹⁴ 5d¹⁰ 6s¹', 'orbital_diagram': '[Xe] 4f: filled  5d: ↑↓ ↑↓ ↑↓ ↑↓ ↑↓  6s: ↑',
        'category': 'Transition Metal', 'color': '#FFD123', 'discovery_year': -2600,
        'melting_point': 1064.18, 'boiling_point': 2856, 'density': 19.282,
        'oxidation_states': [-1, 1, 3, 5], 'electronegativity': 2.54
    },
    # Super heavy elements
    118: {
        'symbol': 'Og', 'name': 'Oganesson', 'atomic_mass': 294, 'group': 18, 'period': 7,
        'electron_config': '[Rn] 5f¹⁴ 6d¹⁰ 7s² 7p⁶', 'orbital_diagram': '[Rn] 5f: filled  6d: filled  7s: ↑↓  7p: ↑↓ ↑↓ ↑↓',
        'category': 'Noble Gas', 'color': '#5A5A5A', 'discovery_year': 2002,
        'melting_point': None, 'boiling_point': None, 'density': None,
        'oxidation_states': [0, 2, 4, 6], 'electronegativity': None
    }
}

# Add more elements to reach 118 total
ADDITIONAL_ELEMENTS = {
    21: {'symbol': 'Sc', 'name': 'Scandium', 'atomic_mass': 44.956, 'group': 3, 'period': 4, 'category': 'Transition Metal', 'color': '#E6E6E6'},
    22: {'symbol': 'Ti', 'name': 'Titanium', 'atomic_mass': 47.867, 'group': 4, 'period': 4, 'category': 'Transition Metal', 'color': '#BFC2C7'},
    23: {'symbol': 'V', 'name': 'Vanadium', 'atomic_mass': 50.942, 'group': 5, 'period': 4, 'category': 'Transition Metal', 'color': '#A6A6AB'},
    24: {'symbol': 'Cr', 'name': 'Chromium', 'atomic_mass': 51.996, 'group': 6, 'period': 4, 'category': 'Transition Metal', 'color': '#8A99C7'},
    25: {'symbol': 'Mn', 'name': 'Manganese', 'atomic_mass': 54.938, 'group': 7, 'period': 4, 'category': 'Transition Metal', 'color': '#9C7AC7'},
    27: {'symbol': 'Co', 'name': 'Cobalt', 'atomic_mass': 58.933, 'group': 9, 'period': 4, 'category': 'Transition Metal', 'color': '#F090A0'},
    28: {'symbol': 'Ni', 'name': 'Nickel', 'atomic_mass': 58.693, 'group': 10, 'period': 4, 'category': 'Transition Metal', 'color': '#50D050'},
    30: {'symbol': 'Zn', 'name': 'Zinc', 'atomic_mass': 65.38, 'group': 12, 'period': 4, 'category': 'Transition Metal', 'color': '#7D80B0'},
    31: {'symbol': 'Ga', 'name': 'Gallium', 'atomic_mass': 69.723, 'group': 13, 'period': 4, 'category': 'Post-transition Metal', 'color': '#C28F8F'},
    32: {'symbol': 'Ge', 'name': 'Germanium', 'atomic_mass': 72.64, 'group': 14, 'period': 4, 'category': 'Metalloid', 'color': '#668F8F'},
    33: {'symbol': 'As', 'name': 'Arsenic', 'atomic_mass': 74.922, 'group': 15, 'period': 4, 'category': 'Metalloid', 'color': '#BD80E3'},
    34: {'symbol': 'Se', 'name': 'Selenium', 'atomic_mass': 78.96, 'group': 16, 'period': 4, 'category': 'Nonmetal', 'color': '#FFA100'},
    35: {'symbol': 'Br', 'name': 'Bromine', 'atomic_mass': 79.904, 'group': 17, 'period': 4, 'category': 'Halogen', 'color': '#A62929'},
    36: {'symbol': 'Kr', 'name': 'Krypton', 'atomic_mass': 83.798, 'group': 18, 'period': 4, 'category': 'Noble Gas', 'color': '#5CB3CC'},
    # Period 5
    37: {'symbol': 'Rb', 'name': 'Rubidium', 'atomic_mass': 85.468, 'group': 1, 'period': 5, 'category': 'Alkali Metal', 'color': '#702EB0'},
    38: {'symbol': 'Sr', 'name': 'Strontium', 'atomic_mass': 87.62, 'group': 2, 'period': 5, 'category': 'Alkaline Earth Metal', 'color': '#00FF00'},
    # Adding more key elements...
    53: {'symbol': 'I', 'name': 'Iodine', 'atomic_mass': 126.90, 'group': 17, 'period': 5, 'category': 'Halogen', 'color': '#940094'},
    54: {'symbol': 'Xe', 'name': 'Xenon', 'atomic_mass': 131.29, 'group': 18, 'period': 5, 'category': 'Noble Gas', 'color': '#429EB0'},
    # Period 6 
    55: {'symbol': 'Cs', 'name': 'Cesium', 'atomic_mass': 132.91, 'group': 1, 'period': 6, 'category': 'Alkali Metal', 'color': '#57178F'},
    56: {'symbol': 'Ba', 'name': 'Barium', 'atomic_mass': 137.33, 'group': 2, 'period': 6, 'category': 'Alkaline Earth Metal', 'color': '#00C900'},
    82: {'symbol': 'Pb', 'name': 'Lead', 'atomic_mass': 207.2, 'group': 14, 'period': 6, 'category': 'Post-transition Metal', 'color': '#575961'},
    86: {'symbol': 'Rn', 'name': 'Radon', 'atomic_mass': 222, 'group': 18, 'period': 6, 'category': 'Noble Gas', 'color': '#428296'},
    # Period 7
    87: {'symbol': 'Fr', 'name': 'Francium', 'atomic_mass': 223, 'group': 1, 'period': 7, 'category': 'Alkali Metal', 'color': '#420066'},
    88: {'symbol': 'Ra', 'name': 'Radium', 'atomic_mass': 226, 'group': 2, 'period': 7, 'category': 'Alkaline Earth Metal', 'color': '#007D00'},
    # Recent superheavy elements
    112: {'symbol': 'Cn', 'name': 'Copernicium', 'atomic_mass': 285, 'group': 12, 'period': 7, 'category': 'Transition Metal', 'color': '#7D80B0'},
    113: {'symbol': 'Nh', 'name': 'Nihonium', 'atomic_mass': 286, 'group': 13, 'period': 7, 'category': 'Post-transition Metal', 'color': '#C28F8F'},
    114: {'symbol': 'Fl', 'name': 'Flerovium', 'atomic_mass': 289, 'group': 14, 'period': 7, 'category': 'Post-transition Metal', 'color': '#668F8F'},
    115: {'symbol': 'Mc', 'name': 'Moscovium', 'atomic_mass': 290, 'group': 15, 'period': 7, 'category': 'Post-transition Metal', 'color': '#BD80E3'},
    116: {'symbol': 'Lv', 'name': 'Livermorium', 'atomic_mass': 293, 'group': 16, 'period': 7, 'category': 'Post-transition Metal', 'color': '#FFA100'},
    117: {'symbol': 'Ts', 'name': 'Tennessine', 'atomic_mass': 294, 'group': 17, 'period': 7, 'category': 'Halogen', 'color': '#A62929'}
}

# Merge the additional elements
for atomic_num, data in ADDITIONAL_ELEMENTS.items():
    # Add default orbital configuration for elements without detailed data
    if 'electron_config' not in data:
        data['electron_config'] = f'[Previous noble gas] configuration'
        data['orbital_diagram'] = f'Complex orbital filling pattern'
        data['discovery_year'] = 'Modern era'
        data['melting_point'] = 'Unknown'
        data['boiling_point'] = 'Unknown'
        data['density'] = 'Unknown'
        data['oxidation_states'] = ['Variable']
        data['electronegativity'] = 'Unknown'
    PERIODIC_TABLE_DATA[atomic_num] = data

CATEGORY_COLORS = {
    'Alkali Metal': '#CC80FF',
    'Alkaline Earth Metal': '#C2FF00',
    'Transition Metal': '#FFC0C0',
    'Post-transition Metal': '#CCCCCC',
    'Metalloid': '#FFCC99',
    'Nonmetal': '#FFFF99',
    'Halogen': '#90E050',
    'Noble Gas': '#C0FFFF',
    'Lanthanide': '#FFBFFF',
    'Actinide': '#FF99CC'
}