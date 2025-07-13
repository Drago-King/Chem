"""
Comprehensive database of commonly used chemical compounds
Organized by categories with IUPAC names, common names, formulas, and SMILES notation.
"""

COMMON_COMPOUNDS_DATABASE = {
    # Organic Acids
    'acetic_acid': {
        'iupac_name': 'ethanoic acid',
        'common_name': 'acetic acid',
        'formula': 'CH₃COOH',
        'smiles': 'CC(=O)O',
        'category': 'Organic Acid',
        'uses': 'Vinegar, food preservative, chemical synthesis',
        'molar_mass': 60.05
    },
    'formic_acid': {
        'iupac_name': 'methanoic acid',
        'common_name': 'formic acid',
        'formula': 'HCOOH',
        'smiles': 'C(=O)O',
        'category': 'Organic Acid',
        'uses': 'Preservative, antibacterial agent',
        'molar_mass': 46.03
    },
    'citric_acid': {
        'iupac_name': '2-hydroxypropane-1,2,3-tricarboxylic acid',
        'common_name': 'citric acid',
        'formula': 'C₆H₈O₇',
        'smiles': 'C(C(=O)O)C(CC(=O)O)(C(=O)O)O',
        'category': 'Organic Acid',
        'uses': 'Food additive, cleaning agent, antioxidant',
        'molar_mass': 192.12
    },
    'lactic_acid': {
        'iupac_name': '2-hydroxypropanoic acid',
        'common_name': 'lactic acid',
        'formula': 'C₃H₆O₃',
        'smiles': 'CC(C(=O)O)O',
        'category': 'Organic Acid',
        'uses': 'Food preservation, cosmetics, biodegradable plastics',
        'molar_mass': 90.08
    },
    'oxalic_acid': {
        'iupac_name': 'ethanedioic acid',
        'common_name': 'oxalic acid',
        'formula': 'C₂H₂O₄',
        'smiles': 'C(=O)(C(=O)O)O',
        'category': 'Organic Acid',
        'uses': 'Rust removal, bleaching agent',
        'molar_mass': 90.03
    },

    # Alcohols
    'methanol': {
        'iupac_name': 'methanol',
        'common_name': 'methyl alcohol',
        'formula': 'CH₃OH',
        'smiles': 'CO',
        'category': 'Alcohol',
        'uses': 'Fuel, solvent, antifreeze',
        'molar_mass': 32.04
    },
    'ethanol': {
        'iupac_name': 'ethanol',
        'common_name': 'ethyl alcohol',
        'formula': 'C₂H₅OH',
        'smiles': 'CCO',
        'category': 'Alcohol',
        'uses': 'Alcoholic beverages, fuel, solvent',
        'molar_mass': 45.08
    },
    'isopropanol': {
        'iupac_name': 'propan-2-ol',
        'common_name': 'isopropyl alcohol',
        'formula': 'C₃H₈O',
        'smiles': 'CC(C)O',
        'category': 'Alcohol',
        'uses': 'Disinfectant, solvent, cleaning agent',
        'molar_mass': 60.10
    },
    'glycerol': {
        'iupac_name': 'propane-1,2,3-triol',
        'common_name': 'glycerin',
        'formula': 'C₃H₈O₃',
        'smiles': 'C(C(CO)O)O',
        'category': 'Alcohol',
        'uses': 'Cosmetics, food additive, antifreeze',
        'molar_mass': 92.09
    },
    'ethylene_glycol': {
        'iupac_name': 'ethane-1,2-diol',
        'common_name': 'ethylene glycol',
        'formula': 'C₂H₆O₂',
        'smiles': 'C(CO)O',
        'category': 'Alcohol',
        'uses': 'Antifreeze, coolant, solvent',
        'molar_mass': 62.07
    },

    # Pharmaceuticals
    'aspirin': {
        'iupac_name': '2-acetoxybenzoic acid',
        'common_name': 'aspirin',
        'formula': 'C₉H₈O₄',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'category': 'Pharmaceutical',
        'uses': 'Pain reliever, anti-inflammatory, blood thinner',
        'molar_mass': 180.16
    },
    'acetaminophen': {
        'iupac_name': 'N-(4-hydroxyphenyl)acetamide',
        'common_name': 'paracetamol',
        'formula': 'C₈H₉NO₂',
        'smiles': 'CC(=O)NC1=CC=C(C=C1)O',
        'category': 'Pharmaceutical',
        'uses': 'Pain reliever, fever reducer',
        'molar_mass': 151.16
    },
    'ibuprofen': {
        'iupac_name': '2-(4-isobutylphenyl)propanoic acid',
        'common_name': 'ibuprofen',
        'formula': 'C₁₃H₁₈O₂',
        'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
        'category': 'Pharmaceutical',
        'uses': 'Anti-inflammatory, pain reliever',
        'molar_mass': 206.28
    },
    'caffeine': {
        'iupac_name': '1,3,7-trimethylpurine-2,6-dione',
        'common_name': 'caffeine',
        'formula': 'C₈H₁₀N₄O₂',
        'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'category': 'Pharmaceutical',
        'uses': 'Stimulant, bronchodilator',
        'molar_mass': 194.19
    },
    'nicotine': {
        'iupac_name': '3-(1-methylpyrrolidin-2-yl)pyridine',
        'common_name': 'nicotine',
        'formula': 'C₁₀H₁₄N₂',
        'smiles': 'CN1CCCC1C2=CN=CC=C2',
        'category': 'Pharmaceutical',
        'uses': 'Stimulant, insecticide',
        'molar_mass': 162.23
    },

    # Inorganic Compounds
    'sodium_chloride': {
        'iupac_name': 'sodium chloride',
        'common_name': 'table salt',
        'formula': 'NaCl',
        'smiles': '[Na+].[Cl-]',
        'category': 'Inorganic Salt',
        'uses': 'Food seasoning, de-icing, chemical production',
        'molar_mass': 58.44
    },
    'calcium_carbonate': {
        'iupac_name': 'calcium carbonate',
        'common_name': 'limestone',
        'formula': 'CaCO₃',
        'smiles': 'C(=O)([O-])[O-].[Ca+2]',
        'category': 'Inorganic Salt',
        'uses': 'Construction, paper making, antacid',
        'molar_mass': 100.09
    },
    'sodium_bicarbonate': {
        'iupac_name': 'sodium hydrogen carbonate',
        'common_name': 'baking soda',
        'formula': 'NaHCO₃',
        'smiles': 'C(=O)(O)[O-].[Na+]',
        'category': 'Inorganic Salt',
        'uses': 'Baking, antacid, fire extinguisher',
        'molar_mass': 84.01
    },
    'potassium_permanganate': {
        'iupac_name': 'potassium permanganate',
        'common_name': 'potassium permanganate',
        'formula': 'KMnO₄',
        'smiles': '[K+].[O-][Mn](=O)(=O)=O',
        'category': 'Inorganic Salt',
        'uses': 'Disinfectant, oxidizing agent',
        'molar_mass': 158.03
    },
    'ammonium_nitrate': {
        'iupac_name': 'ammonium nitrate',
        'common_name': 'ammonium nitrate',
        'formula': 'NH₄NO₃',
        'smiles': '[NH4+].[O-][N+](=O)[O-]',
        'category': 'Inorganic Salt',
        'uses': 'Fertilizer, explosive',
        'molar_mass': 80.04
    },

    # Vitamins
    'ascorbic_acid': {
        'iupac_name': '(5R)-5-[(1S)-1,2-dihydroxyethyl]-3,4-dihydroxyfuran-2(5H)-one',
        'common_name': 'vitamin C',
        'formula': 'C₆H₈O₆',
        'smiles': 'C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O',
        'category': 'Vitamin',
        'uses': 'Antioxidant, immune system support',
        'molar_mass': 176.12
    },
    'thiamine': {
        'iupac_name': '2-[3-[(4-amino-2-methylpyrimidin-5-yl)methyl]-4-methyl-1,3-thiazol-3-ium-5-yl]ethanol',
        'common_name': 'vitamin B1',
        'formula': 'C₁₂H₁₇N₄OS⁺',
        'smiles': 'CC1=NC(=C(N=C1N)C)[N+]2=CSC(=C2)CCO',
        'category': 'Vitamin',
        'uses': 'Metabolism, nervous system function',
        'molar_mass': 265.35
    },
    'retinol': {
        'iupac_name': '(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohex-1-en-1-yl)nona-2,4,6,8-tetraen-1-ol',
        'common_name': 'vitamin A',
        'formula': 'C₂₀H₃₀O',
        'smiles': 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C',
        'category': 'Vitamin',
        'uses': 'Vision, immune function, cell growth',
        'molar_mass': 286.45
    },

    # Polymers and Plastics
    'polyethylene': {
        'iupac_name': 'poly(ethene)',
        'common_name': 'polyethylene',
        'formula': '(C₂H₄)ₙ',
        'smiles': 'CCCCCCCC',
        'category': 'Polymer',
        'uses': 'Plastic bags, bottles, packaging',
        'molar_mass': 28.05
    },
    'polystyrene': {
        'iupac_name': 'poly(1-phenylethene)',
        'common_name': 'polystyrene',
        'formula': '(C₈H₈)ₙ',
        'smiles': 'CC(C1=CC=CC=C1)C(C2=CC=CC=C2)C',
        'category': 'Polymer',
        'uses': 'Disposable cups, insulation, packaging',
        'molar_mass': 104.15
    },
    'nylon': {
        'iupac_name': 'poly(hexamethylene adipamide)',
        'common_name': 'nylon-6,6',
        'formula': '(C₁₂H₂₂N₂O₂)ₙ',
        'smiles': 'C(CCCCC)NCCCCCCNC(=O)CCCCC(=O)',
        'category': 'Polymer',
        'uses': 'Textiles, rope, mechanical parts',
        'molar_mass': 226.32
    },

    # Food Compounds
    'glucose': {
        'iupac_name': '(2R,3S,4R,5R)-2,3,4,5,6-pentahydroxyhexanal',
        'common_name': 'glucose',
        'formula': 'C₆H₁₂O₆',
        'smiles': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O',
        'category': 'Sugar',
        'uses': 'Energy source, sweetener',
        'molar_mass': 180.16
    },
    'fructose': {
        'iupac_name': '(3S,4R,5R)-1,3,4,5,6-pentahydroxyhexan-2-one',
        'common_name': 'fructose',
        'formula': 'C₆H₁₂O₆',
        'smiles': 'C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O',
        'category': 'Sugar',
        'uses': 'Sweetener, energy source',
        'molar_mass': 180.16
    },
    'sucrose': {
        'iupac_name': '(2R,3R,4S,5S,6R)-2-[(2S,3S,4S,5R)-3,4-dihydroxy-2,5-bis(hydroxymethyl)oxolan-2-yl]oxy-6-(hydroxymethyl)oxane-3,4,5-triol',
        'common_name': 'table sugar',
        'formula': 'C₁₂H₂₂O₁₁',
        'smiles': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
        'category': 'Sugar',
        'uses': 'Sweetener, food preservation',
        'molar_mass': 342.30
    },
    'lactose': {
        'iupac_name': '(2R,3R,4S,5R,6S)-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6S)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-oxane-3,4,5-triol',
        'common_name': 'milk sugar',
        'formula': 'C₁₂H₂₂O₁₁',
        'smiles': 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O',
        'category': 'Sugar',
        'uses': 'Milk component, pharmaceutical excipient',
        'molar_mass': 342.30
    },

    # Industrial Chemicals
    'sulfuric_acid': {
        'iupac_name': 'sulfuric acid',
        'common_name': 'oil of vitriol',
        'formula': 'H₂SO₄',
        'smiles': 'OS(=O)(=O)O',
        'category': 'Inorganic Acid',
        'uses': 'Chemical manufacturing, batteries, fertilizers',
        'molar_mass': 98.08
    },
    'hydrochloric_acid': {
        'iupac_name': 'hydrogen chloride',
        'common_name': 'muriatic acid',
        'formula': 'HCl',
        'smiles': 'Cl',
        'category': 'Inorganic Acid',
        'uses': 'Metal cleaning, pH control, food processing',
        'molar_mass': 36.46
    },
    'nitric_acid': {
        'iupac_name': 'nitric acid',
        'common_name': 'aqua fortis',
        'formula': 'HNO₃',
        'smiles': '[N+](=O)(O)[O-]',
        'category': 'Inorganic Acid',
        'uses': 'Fertilizers, explosives, metal etching',
        'molar_mass': 63.01
    },
    'ammonia': {
        'iupac_name': 'ammonia',
        'common_name': 'ammonia',
        'formula': 'NH₃',
        'smiles': 'N',
        'category': 'Inorganic Compound',
        'uses': 'Fertilizers, cleaning agents, refrigerant',
        'molar_mass': 17.03
    },
    'hydrogen_peroxide': {
        'iupac_name': 'hydrogen peroxide',
        'common_name': 'hydrogen peroxide',
        'formula': 'H₂O₂',
        'smiles': 'OO',
        'category': 'Inorganic Compound',
        'uses': 'Disinfectant, bleaching agent, rocket fuel',
        'molar_mass': 34.01
    },

    # Organic Solvents
    'acetone': {
        'iupac_name': 'propan-2-one',
        'common_name': 'acetone',
        'formula': 'C₃H₆O',
        'smiles': 'CC(=O)C',
        'category': 'Ketone',
        'uses': 'Solvent, nail polish remover, plastic production',
        'molar_mass': 58.08
    },
    'benzene': {
        'iupac_name': 'benzene',
        'common_name': 'benzene',
        'formula': 'C₆H₆',
        'smiles': 'C1=CC=CC=C1',
        'category': 'Aromatic Hydrocarbon',
        'uses': 'Chemical synthesis, gasoline additive',
        'molar_mass': 78.11
    },
    'toluene': {
        'iupac_name': 'methylbenzene',
        'common_name': 'toluene',
        'formula': 'C₇H₈',
        'smiles': 'CC1=CC=CC=C1',
        'category': 'Aromatic Hydrocarbon',
        'uses': 'Solvent, fuel additive, chemical synthesis',
        'molar_mass': 92.14
    },
    'dichloromethane': {
        'iupac_name': 'dichloromethane',
        'common_name': 'methylene chloride',
        'formula': 'CH₂Cl₂',
        'smiles': 'C(Cl)Cl',
        'category': 'Halogenated Solvent',
        'uses': 'Paint stripper, degreaser, extraction solvent',
        'molar_mass': 84.93
    },
    'chloroform': {
        'iupac_name': 'trichloromethane',
        'common_name': 'chloroform',
        'formula': 'CHCl₃',
        'smiles': 'C(Cl)(Cl)Cl',
        'category': 'Halogenated Solvent',
        'uses': 'Solvent, refrigerant precursor',
        'molar_mass': 119.38
    },

    # Biochemical Compounds
    'cholesterol': {
        'iupac_name': '(3S,8S,9S,10R,13R,14S,17R)-10,13-dimethyl-17-[(2R)-6-methylheptan-2-yl]-2,3,4,7,8,9,11,12,14,15,16,17-dodecahydro-1H-cyclopenta[a]phenanthren-3-ol',
        'common_name': 'cholesterol',
        'formula': 'C₂₇H₄₆O',
        'smiles': 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
        'category': 'Steroid',
        'uses': 'Cell membrane component, hormone precursor',
        'molar_mass': 386.65
    },
    'alanine': {
        'iupac_name': '2-aminopropanoic acid',
        'common_name': 'alanine',
        'formula': 'C₃H₇NO₂',
        'smiles': 'C[C@@H](C(=O)O)N',
        'category': 'Amino Acid',
        'uses': 'Protein synthesis, energy source',
        'molar_mass': 89.09
    },
    'glycine': {
        'iupac_name': '2-aminoacetic acid',
        'common_name': 'glycine',
        'formula': 'C₂H₅NO₂',
        'smiles': 'C(C(=O)O)N',
        'category': 'Amino Acid',
        'uses': 'Protein synthesis, neurotransmitter',
        'molar_mass': 75.07
    },
    'urea': {
        'iupac_name': 'urea',
        'common_name': 'carbamide',
        'formula': 'CH₄N₂O',
        'smiles': 'C(=O)(N)N',
        'category': 'Organic Compound',
        'uses': 'Fertilizer, cosmetics, animal feed',
        'molar_mass': 60.06
    },

    # Gases
    'carbon_dioxide': {
        'iupac_name': 'carbon dioxide',
        'common_name': 'carbon dioxide',
        'formula': 'CO₂',
        'smiles': 'C(=O)=O',
        'category': 'Inorganic Gas',
        'uses': 'Fire extinguisher, dry ice, carbonation',
        'molar_mass': 44.01
    },
    'carbon_monoxide': {
        'iupac_name': 'carbon monoxide',
        'common_name': 'carbon monoxide',
        'formula': 'CO',
        'smiles': '[C-]#[O+]',
        'category': 'Inorganic Gas',
        'uses': 'Reducing agent, fuel gas',
        'molar_mass': 28.01
    },
    'methane': {
        'iupac_name': 'methane',
        'common_name': 'natural gas',
        'formula': 'CH₄',
        'smiles': 'C',
        'category': 'Alkane',
        'uses': 'Fuel, heating, chemical feedstock',
        'molar_mass': 16.04
    },
    'ethane': {
        'iupac_name': 'ethane',
        'common_name': 'ethane',
        'formula': 'C₂H₆',
        'smiles': 'CC',
        'category': 'Alkane',
        'uses': 'Fuel, petrochemical feedstock',
        'molar_mass': 30.07
    },
    'propane': {
        'iupac_name': 'propane',
        'common_name': 'propane',
        'formula': 'C₃H₈',
        'smiles': 'CCC',
        'category': 'Alkane',
        'uses': 'Fuel, heating, refrigerant',
        'molar_mass': 44.10
    },
    'butane': {
        'iupac_name': 'butane',
        'common_name': 'butane',
        'formula': 'C₄H₁₀',
        'smiles': 'CCCC',
        'category': 'Alkane',
        'uses': 'Fuel, lighter fluid, aerosol propellant',
        'molar_mass': 58.12
    },

    # Additional Important Compounds
    'sodium_hydroxide': {
        'iupac_name': 'sodium hydroxide',
        'common_name': 'caustic soda',
        'formula': 'NaOH',
        'smiles': '[OH-].[Na+]',
        'category': 'Inorganic Base',
        'uses': 'Soap making, drain cleaner, chemical manufacturing',
        'molar_mass': 40.00
    },
    'calcium_oxide': {
        'iupac_name': 'calcium oxide',
        'common_name': 'quicklime',
        'formula': 'CaO',
        'smiles': '[Ca+2].[O-2]',
        'category': 'Inorganic Oxide',
        'uses': 'Cement, steel production, water treatment',
        'molar_mass': 56.08
    },
    'silicon_dioxide': {
        'iupac_name': 'silicon dioxide',
        'common_name': 'silica',
        'formula': 'SiO₂',
        'smiles': 'O=[Si]=O',
        'category': 'Inorganic Oxide',
        'uses': 'Glass making, electronics, construction',
        'molar_mass': 60.08
    },
    'aluminum_oxide': {
        'iupac_name': 'aluminum oxide',
        'common_name': 'alumina',
        'formula': 'Al₂O₃',
        'smiles': '[O-2].[O-2].[O-2].[Al+3].[Al+3]',
        'category': 'Inorganic Oxide',
        'uses': 'Abrasives, ceramics, aluminum production',
        'molar_mass': 101.96
    },
    'iron_oxide': {
        'iupac_name': 'iron(III) oxide',
        'common_name': 'rust',
        'formula': 'Fe₂O₃',
        'smiles': '[O-2].[O-2].[O-2].[Fe+3].[Fe+3]',
        'category': 'Inorganic Oxide',
        'uses': 'Pigments, polishing compounds, thermite',
        'molar_mass': 159.69
    },

    # Additional Organic Compounds
    'benzaldehyde': {
        'iupac_name': 'benzaldehyde',
        'common_name': 'benzaldehyde',
        'formula': 'C₇H₆O',
        'smiles': 'C1=CC=C(C=C1)C=O',
        'category': 'Aldehyde',
        'uses': 'Flavoring agent, perfumes, dyes',
        'molar_mass': 106.12
    },
    'formaldehyde': {
        'iupac_name': 'methanal',
        'common_name': 'formaldehyde',
        'formula': 'CH₂O',
        'smiles': 'C=O',
        'category': 'Aldehyde',
        'uses': 'Disinfectant, preservative, resins',
        'molar_mass': 30.03
    },
    'acetaldehyde': {
        'iupac_name': 'ethanal',
        'common_name': 'acetaldehyde',
        'formula': 'C₂H₄O',
        'smiles': 'CC=O',
        'category': 'Aldehyde',
        'uses': 'Chemical synthesis, flavoring',
        'molar_mass': 44.05
    },

    # More Pharmaceuticals
    'morphine': {
        'iupac_name': '(5α,6α)-7,8-didehydro-4,5-epoxy-17-methylmorphinan-3,6-diol',
        'common_name': 'morphine',
        'formula': 'C₁₇H₁₉NO₃',
        'smiles': 'CN1CC[C@]23C4=C5C(=C(C=C4[C@H]1[C@@H]2C=C[C@@H]3O)O)OC=C5',
        'category': 'Pharmaceutical',
        'uses': 'Pain management, opioid analgesic',
        'molar_mass': 285.34
    },
    'penicillin': {
        'iupac_name': '(2S,5R,6R)-3,3-dimethyl-7-oxo-6-[(phenylacetyl)amino]-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid',
        'common_name': 'penicillin G',
        'formula': 'C₁₆H₁₈N₂O₄S',
        'smiles': 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
        'category': 'Pharmaceutical',
        'uses': 'Antibiotic, bacterial infections',
        'molar_mass': 334.39
    },
    'codeine': {
        'iupac_name': '(5α,6α)-7,8-didehydro-4,5-epoxy-3-methoxy-17-methylmorphinan-6-ol',
        'common_name': 'codeine',
        'formula': 'C₁₈H₂₁NO₃',
        'smiles': 'COC1=C(C=C2C[C@H]3N(CC[C@@]24C1=C(C=C4)[C@H]([C@@H]3O)O)C)O',
        'category': 'Pharmaceutical',
        'uses': 'Pain relief, cough suppressant',
        'molar_mass': 299.36
    },
    'insulin': {
        'iupac_name': 'insulin',
        'common_name': 'insulin',
        'formula': 'C₂₅₇H₃₈₃N₆₅O₇₇S₆',
        'smiles': 'CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC2=CNC3=CC=CC=C32)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC4=CC=C(C=C4)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC5=CC=CC=C5)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]6CCCN6C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC7=CC=CC=C7)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC8=CC=C(C=C8)O)NC(=O)[C@@H]9CCCN9C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%10=CC=CC=C%10)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%11=CC=C(C=C%11)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%12=CNC%13=CC=CC=C%12%13)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%14=CC=CC=C%14)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%15CCCN%15C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%16=CC=CC=C%16)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%17=CC=C(C=C%17)O)NC(=O)[C@@H]%18CCCN%18C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%19=CC=CC=C%19)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%20=CC=C(C=C%20)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%21=CNC%22=CC=CC=C%21%22)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%23=CC=CC=C%23)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%24CCCN%24C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%25=CC=CC=C%25)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%26=CC=C(C=C%26)O)NC(=O)[C@@H]%27CCCN%27C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%28=CC=CC=C%28)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%29=CC=C(C=C%29)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%30=CNC%31=CC=CC=C%30%31)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%32=CC=CC=C%32)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%33CCCN%33C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%34=CC=CC=C%34)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%35=CC=C(C=C%35)O)NC(=O)[C@@H]%36CCCN%36C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%37=CC=CC=C%37)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%38=CC=C(C=C%38)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%39=CNC%40=CC=CC=C%39%40)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%41=CC=CC=C%41)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%42CCCN%42C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%43=CC=CC=C%43)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%44=CC=C(C=C%44)O)NC(=O)[C@@H]%45CCCN%45C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%46=CC=CC=C%46)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%47=CC=C(C=C%47)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%48=CNC%49=CC=CC=C%48%49)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%50=CC=CC=C%50)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%51CCCN%51C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%52=CC=CC=C%52)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%53=CC=C(C=C%53)O)NC(=O)[C@@H]%54CCCN%54C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%55=CC=CC=C%55)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%56=CC=C(C=C%56)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%57=CNC%58=CC=CC=C%57%58)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%59=CC=CC=C%59)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%60CCCN%60C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%61=CC=CC=C%61)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%62=CC=C(C=C%62)O)NC(=O)[C@@H]%63CCCN%63C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%64=CC=CC=C%64)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%65=CC=C(C=C%65)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%66=CNC%67=CC=CC=C%66%67)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%68=CC=CC=C%68)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%69CCCN%69C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%70=CC=CC=C%70)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%71=CC=C(C=C%71)O)NC(=O)[C@@H]%72CCCN%72C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%73=CC=CC=C%73)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%74=CC=C(C=C%74)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%75=CNC%76=CC=CC=C%75%76)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%77=CC=CC=C%77)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%78CCCN%78C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%79=CC=CC=C%79)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%80=CC=C(C=C%80)O)NC(=O)[C@@H]%81CCCN%81C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%82=CC=CC=C%82)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%83=CC=C(C=C%83)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%84=CNC%85=CC=CC=C%84%85)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%86=CC=CC=C%86)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%87CCCN%87C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%88=CC=CC=C%88)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%89=CC=C(C=C%89)O)NC(=O)[C@@H]%90CCCN%90C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%91=CC=CC=C%91)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%92=CC=C(C=C%92)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%93=CNC%94=CC=CC=C%93%94)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%95=CC=CC=C%95)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%96CCCN%96C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%97=CC=CC=C%97)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%98=CC=C(C=C%98)O)NC(=O)[C@@H]%99CCCN%99C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%100=CC=CC=C%100)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%101=CC=C(C=C%101)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%102=CNC%103=CC=CC=C%102%103)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%104=CC=CC=C%104)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%105CCCN%105C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%106=CC=CC=C%106)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%107=CC=C(C=C%107)O)NC(=O)[C@@H]%108CCCN%108C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%109=CC=CC=C%109)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%110=CC=C(C=C%110)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%111=CNC%112=CC=CC=C%111%112)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%113=CC=CC=C%113)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%114CCCN%114C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%115=CC=CC=C%115)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%116=CC=C(C=C%116)O)NC(=O)[C@@H]%117CCCN%117C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%118=CC=CC=C%118)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%119=CC=C(C=C%119)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%120=CNC%121=CC=CC=C%120%121)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%122=CC=CC=C%122)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%123CCCN%123C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%124=CC=CC=C%124)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%125=CC=C(C=C%125)O)NC(=O)[C@@H]%126CCCN%126C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%127=CC=CC=C%127)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%128=CC=C(C=C%128)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%129=CNC%130=CC=CC=C%129%130)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%131=CC=CC=C%131)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%132CCCN%132C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%133=CC=CC=C%133)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%134=CC=C(C=C%134)O)NC(=O)[C@@H]%135CCCN%135C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%136=CC=CC=C%136)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%137=CC=C(C=C%137)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%138=CNC%139=CC=CC=C%138%139)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%140=CC=CC=C%140)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%141CCCN%141C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%142=CC=CC=C%142)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%143=CC=C(C=C%143)O)NC(=O)[C@@H]%144CCCN%144C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%145=CC=CC=C%145)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%146=CC=C(C=C%146)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%147=CNC%148=CC=CC=C%147%148)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%149=CC=CC=C%149)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%150CCCN%150C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%151=CC=CC=C%151)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%152=CC=C(C=C%152)O)NC(=O)[C@@H]%153CCCN%153C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%154=CC=CC=C%154)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%155=CC=C(C=C%155)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%156=CNC%157=CC=CC=C%156%157)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%158=CC=CC=C%158)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%159CCCN%159C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%160=CC=CC=C%160)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%161=CC=C(C=C%161)O)NC(=O)[C@@H]%162CCCN%162C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%163=CC=CC=C%163)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%164=CC=C(C=C%164)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%165=CNC%166=CC=CC=C%165%166)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%167=CC=CC=C%167)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%168CCCN%168C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%169=CC=CC=C%169)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%170=CC=C(C=C%170)O)NC(=O)[C@@H]%171CCCN%171C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%172=CC=CC=C%172)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%173=CC=C(C=C%173)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%174=CNC%175=CC=CC=C%174%175)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%176=CC=CC=C%176)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%177CCCN%177C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%178=CC=CC=C%178)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%179=CC=C(C=C%179)O)NC(=O)[C@@H]%180CCCN%180C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%181=CC=CC=C%181)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%182=CC=C(C=C%182)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%183=CNC%184=CC=CC=C%183%184)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%185=CC=CC=C%185)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%186CCCN%186C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%187=CC=CC=C%187)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%188=CC=C(C=C%188)O)NC(=O)[C@@H]%189CCCN%189C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%190=CC=CC=C%190)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%191=CC=C(C=C%191)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%192=CNC%193=CC=CC=C%192%193)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%194=CC=CC=C%194)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%195CCCN%195C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%196=CC=CC=C%196)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%197=CC=C(C=C%197)O)NC(=O)[C@@H]%198CCCN%198C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%199=CC=CC=C%199)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%200=CC=C(C=C%200)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%201=CNC%202=CC=CC=C%201%202)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%203=CC=CC=C%203)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%204CCCN%204C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%205=CC=CC=C%205)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%206=CC=C(C=C%206)O)NC(=O)[C@@H]%207CCCN%207C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%208=CC=CC=C%208)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%209=CC=C(C=C%209)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%210=CNC%211=CC=CC=C%210%211)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%212=CC=CC=C%212)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%213CCCN%213C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%214=CC=CC=C%214)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%215=CC=C(C=C%215)O)NC(=O)[C@@H]%216CCCN%216C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%217=CC=CC=C%217)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%218=CC=C(C=C%218)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%219=CNC%220=CC=CC=C%219%220)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%221=CC=CC=C%221)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%222CCCN%222C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%223=CC=CC=C%223)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%224=CC=C(C=C%224)O)NC(=O)[C@@H]%225CCCN%225C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%226=CC=CC=C%226)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%227=CC=C(C=C%227)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%228=CNC%229=CC=CC=C%228%229)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%230=CC=CC=C%230)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%231CCCN%231C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%232=CC=CC=C%232)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%233=CC=C(C=C%233)O)NC(=O)[C@@H]%234CCCN%234C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%235=CC=CC=C%235)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%236=CC=C(C=C%236)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%237=CNC%238=CC=CC=C%237%238)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%239=CC=CC=C%239)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%240CCCN%240C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%241=CC=CC=C%241)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%242=CC=C(C=C%242)O)NC(=O)[C@@H]%243CCCN%243C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%244=CC=CC=C%244)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%245=CC=C(C=C%245)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%246=CNC%247=CC=CC=C%246%247)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%248=CC=CC=C%248)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%249CCCN%249C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%250=CC=CC=C%250)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%251=CC=C(C=C%251)O)NC(=O)[C@@H]%252CCCN%252C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%253=CC=CC=C%253)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%254=CC=C(C=C%254)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%255=CNC%256=CC=CC=C%255%256)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%257=CC=CC=C%257)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%258CCCN%258C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%259=CC=CC=C%259)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%260=CC=C(C=C%260)O)NC(=O)[C@@H]%261CCCN%261C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%262=CC=CC=C%262)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%263=CC=C(C=C%263)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%264=CNC%265=CC=CC=C%264%265)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%266=CC=CC=C%266)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%267CCCN%267C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%268=CC=CC=C%268)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%269=CC=C(C=C%269)O)NC(=O)[C@@H]%270CCCN%270C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%271=CC=CC=C%271)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%272=CC=C(C=C%272)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%273=CNC%274=CC=CC=C%273%274)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%275=CC=CC=C%275)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%276CCCN%276C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%277=CC=CC=C%277)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%278=CC=C(C=C%278)O)NC(=O)[C@@H]%279CCCN%279C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%280=CC=CC=C%280)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%281=CC=C(C=C%281)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%282=CNC%283=CC=CC=C%282%283)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%284=CC=CC=C%284)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%285CCCN%285C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%286=CC=CC=C%286)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%287=CC=C(C=C%287)O)NC(=O)[C@@H]%288CCCN%288C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%289=CC=CC=C%289)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%290=CC=C(C=C%290)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%291=CNC%292=CC=CC=C%291%292)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%293=CC=CC=C%293)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%294CCCN%294C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%295=CC=CC=C%295)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%296=CC=C(C=C%296)O)NC(=O)[C@@H]%297CCCN%297C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%298=CC=CC=C%298)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%299=CC=C(C=C%299)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%300=CNC%301=CC=CC=C%300%301)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%302=CC=CC=C%302)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%303CCCN%303C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%304=CC=CC=C%304)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%305=CC=C(C=C%305)O)NC(=O)[C@@H]%306CCCN%306C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%307=CC=CC=C%307)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%308=CC=C(C=C%308)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%309=CNC%310=CC=CC=C%309%310)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%311=CC=CC=C%311)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%312CCCN%312C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%313=CC=CC=C%313)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%314=CC=C(C=C%314)O)NC(=O)[C@@H]%315CCCN%315C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%316=CC=CC=C%316)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%317=CC=C(C=C%317)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%318=CNC%319=CC=CC=C%318%319)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%320=CC=CC=C%320)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%321CCCN%321C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%322=CC=CC=C%322)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%323=CC=C(C=C%323)O)NC(=O)[C@@H]%324CCCN%324C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%325=CC=CC=C%325)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%326=CC=C(C=C%326)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%327=CNC%328=CC=CC=C%327%328)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%329=CC=CC=C%329)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%330CCCN%330C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%331=CC=CC=C%331)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%332=CC=C(C=C%332)O)NC(=O)[C@@H]%333CCCN%333C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%334=CC=CC=C%334)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%335=CC=C(C=C%335)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%336=CNC%337=CC=CC=C%336%337)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%338=CC=CC=C%338)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%339CCCN%339C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%340=CC=CC=C%340)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%341=CC=C(C=C%341)O)NC(=O)[C@@H]%342CCCN%342C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%343=CC=CC=C%343)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%344=CC=C(C=C%344)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%345=CNC%346=CC=CC=C%345%346)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%347=CC=CC=C%347)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%348CCCN%348C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%349=CC=CC=C%349)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%350=CC=C(C=C%350)O)NC(=O)[C@@H]%351CCCN%351C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%352=CC=CC=C%352)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%353=CC=C(C=C%353)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%354=CNC%355=CC=CC=C%354%355)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%356=CC=CC=C%356)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%357CCCN%357C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%358=CC=CC=C%358)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%359=CC=C(C=C%359)O)NC(=O)[C@@H]%360CCCN%360C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%361=CC=CC=C%361)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%362=CC=C(C=C%362)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%363=CNC%364=CC=CC=C%363%364)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%365=CC=CC=C%365)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%366CCCN%366C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%367=CC=CC=C%367)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%368=CC=C(C=C%368)O)NC(=O)[C@@H]%369CCCN%369C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%370=CC=CC=C%370)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%371=CC=C(C=C%371)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%372=CNC%373=CC=CC=C%372%373)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%374=CC=CC=C%374)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]%375CCCN%375C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](C)NC(=O)[C@H](CC%376=CC=CC=C%376)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC%377=CC=C(C=C%377)O)NC(=O)[C@@H]%378CCCN%378C(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%379=CC=CC=C%379)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CC%380=CC=C(C=C%380)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CC%381=CNC%382=CC=CC=C%381%382)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC%383=CC=CC=C%383)C(=O)O',
        'category': 'Pharmaceutical',
        'uses': 'Diabetes treatment, blood sugar regulation',
        'molar_mass': 5808
    },
    'dopamine': {
        'iupac_name': '4-(2-aminoethyl)benzene-1,2-diol',
        'common_name': 'dopamine',
        'formula': 'C₈H₁₁NO₂',
        'smiles': 'C1=CC(=C(C=C1CCN)O)O',
        'category': 'Pharmaceutical',
        'uses': 'Neurotransmitter, Parkinson\'s disease treatment',
        'molar_mass': 153.18
    },

    # More Esters and Ethers
    'ethyl_acetate': {
        'iupac_name': 'ethyl acetate',
        'common_name': 'ethyl acetate',
        'formula': 'C₄H₈O₂',
        'smiles': 'CCOC(=O)C',
        'category': 'Ester',
        'uses': 'Solvent, nail polish remover, glues',
        'molar_mass': 88.11
    },
    'methyl_salicylate': {
        'iupac_name': 'methyl 2-hydroxybenzoate',
        'common_name': 'oil of wintergreen',
        'formula': 'C₈H₈O₃',
        'smiles': 'COC(=O)C1=CC=CC=C1O',
        'category': 'Ester',
        'uses': 'Flavoring, topical analgesic, perfumes',
        'molar_mass': 152.15
    },
    'diethyl_ether': {
        'iupac_name': 'diethyl ether',
        'common_name': 'ether',
        'formula': 'C₄H₁₀O',
        'smiles': 'CCOCC',
        'category': 'Ether',
        'uses': 'Anesthetic, solvent, fuel additive',
        'molar_mass': 74.12
    },

    # More Hydrocarbons
    'xylene': {
        'iupac_name': 'dimethylbenzene',
        'common_name': 'xylene',
        'formula': 'C₈H₁₀',
        'smiles': 'CC1=CC=CC=C1C',
        'category': 'Aromatic Hydrocarbon',
        'uses': 'Solvent, gasoline additive, chemical synthesis',
        'molar_mass': 106.17
    },
    'styrene': {
        'iupac_name': 'phenylethene',
        'common_name': 'styrene',
        'formula': 'C₈H₈',
        'smiles': 'C=CC1=CC=CC=C1',
        'category': 'Aromatic Hydrocarbon',
        'uses': 'Plastic production, synthetic rubber',
        'molar_mass': 104.15
    },
    'naphthalene': {
        'iupac_name': 'naphthalene',
        'common_name': 'naphthalene',
        'formula': 'C₁₀H₈',
        'smiles': 'C1=CC=C2C=CC=CC2=C1',
        'category': 'Aromatic Hydrocarbon',
        'uses': 'Mothballs, chemical synthesis, preservative',
        'molar_mass': 128.17
    },

    # More Industrial Chemicals
    'phosphoric_acid': {
        'iupac_name': 'phosphoric acid',
        'common_name': 'orthophosphoric acid',
        'formula': 'H₃PO₄',
        'smiles': 'OP(=O)(O)O',
        'category': 'Inorganic Acid',
        'uses': 'Food additive, fertilizers, rust removal',
        'molar_mass': 97.99
    },
    'sodium_carbonate': {
        'iupac_name': 'sodium carbonate',
        'common_name': 'soda ash',
        'formula': 'Na₂CO₃',
        'smiles': 'C(=O)([O-])[O-].[Na+].[Na+]',
        'category': 'Inorganic Salt',
        'uses': 'Glass making, detergents, water treatment',
        'molar_mass': 105.99
    },
    'potassium_chloride': {
        'iupac_name': 'potassium chloride',
        'common_name': 'potash',
        'formula': 'KCl',
        'smiles': '[K+].[Cl-]',
        'category': 'Inorganic Salt',
        'uses': 'Fertilizer, salt substitute, medical uses',
        'molar_mass': 74.55
    },

    # More Biomolecules
    'phenylalanine': {
        'iupac_name': '2-amino-3-phenylpropanoic acid',
        'common_name': 'phenylalanine',
        'formula': 'C₉H₁₁NO₂',
        'smiles': 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
        'category': 'Amino Acid',
        'uses': 'Protein synthesis, precursor to neurotransmitters',
        'molar_mass': 165.19
    },
    'leucine': {
        'iupac_name': '2-amino-4-methylpentanoic acid',
        'common_name': 'leucine',
        'formula': 'C₆H₁₃NO₂',
        'smiles': 'CC(C)C[C@@H](C(=O)O)N',
        'category': 'Amino Acid',
        'uses': 'Protein synthesis, muscle building',
        'molar_mass': 131.17
    },
    'lysine': {
        'iupac_name': '2,6-diaminohexanoic acid',
        'common_name': 'lysine',
        'formula': 'C₆H₁₄N₂O₂',
        'smiles': 'C(CCN)C[C@@H](C(=O)O)N',
        'category': 'Amino Acid',
        'uses': 'Protein synthesis, essential amino acid',
        'molar_mass': 146.19
    },

    # More Vitamins and Nutrients
    'riboflavin': {
        'iupac_name': '7,8-dimethyl-10-[(2S,3S,4R)-2,3,4,5-tetrahydroxypentyl]benzo[g]pteridine-2,4-dione',
        'common_name': 'vitamin B2',
        'formula': 'C₁₇H₂₀N₄O₆',
        'smiles': 'CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@H]([C@H]([C@H](CO)O)O)O',
        'category': 'Vitamin',
        'uses': 'Energy metabolism, antioxidant',
        'molar_mass': 376.36
    },
    'niacin': {
        'iupac_name': 'pyridine-3-carboxylic acid',
        'common_name': 'vitamin B3',
        'formula': 'C₆H₅NO₂',
        'smiles': 'C1=CC(=CN=C1)C(=O)O',
        'category': 'Vitamin',
        'uses': 'Cholesterol management, energy metabolism',
        'molar_mass': 123.11
    },
    'folate': {
        'iupac_name': 'N-[4-[[(2-amino-4-oxo-1H-pteridin-6-yl)methyl]amino]benzoyl]-L-glutamic acid',
        'common_name': 'folic acid',
        'formula': 'C₁₉H₁₉N₇O₆',
        'smiles': 'C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N',
        'category': 'Vitamin',
        'uses': 'DNA synthesis, red blood cell formation',
        'molar_mass': 441.40
    },

    # Additional Food Compounds
    'vanillin': {
        'iupac_name': '4-hydroxy-3-methoxybenzaldehyde',
        'common_name': 'vanillin',
        'formula': 'C₈H₈O₃',
        'smiles': 'COC1=C(C=CC(=C1)C=O)O',
        'category': 'Aldehyde',
        'uses': 'Flavoring agent, food additive, perfumes',
        'molar_mass': 152.15
    },
    'menthol': {
        'iupac_name': '(1R,2S,5R)-2-isopropyl-5-methylcyclohexan-1-ol',
        'common_name': 'menthol',
        'formula': 'C₁₀H₂₀O',
        'smiles': 'C[C@H]1CC[C@@H]([C@H](C1)O)C(C)C',
        'category': 'Alcohol',
        'uses': 'Flavoring, cooling agent, topical analgesic',
        'molar_mass': 156.27
    },
    'thymol': {
        'iupac_name': '2-isopropyl-5-methylphenol',
        'common_name': 'thymol',
        'formula': 'C₁₀H₁₄O',
        'smiles': 'CC1=CC=C(C(=C1)C(C)C)O',
        'category': 'Phenol',
        'uses': 'Antiseptic, fungicide, flavoring agent',
        'molar_mass': 150.22
    },

    # More Environmental/Energy Compounds
    'freon_12': {
        'iupac_name': 'dichlorodifluoromethane',
        'common_name': 'freon-12',
        'formula': 'CCl₂F₂',
        'smiles': 'C(Cl)(Cl)(F)F',
        'category': 'Refrigerant',
        'uses': 'Refrigerant (now banned), aerosol propellant',
        'molar_mass': 120.91
    },
    'ddt': {
        'iupac_name': '1,1,1-trichloro-2,2-bis(4-chlorophenyl)ethane',
        'common_name': 'DDT',
        'formula': 'C₁₄H₉Cl₅',
        'smiles': 'C(C1=CC=C(C=C1)Cl)(C2=CC=C(C=C2)Cl)C(Cl)(Cl)Cl',
        'category': 'Pesticide',
        'uses': 'Insecticide (now banned in many countries)',
        'molar_mass': 354.49
    },

    # Additional Essential Compounds to reach 100+
    'pentane': {
        'iupac_name': 'pentane',
        'common_name': 'pentane',
        'formula': 'C₅H₁₂',
        'smiles': 'CCCCC',
        'category': 'Alkane',
        'uses': 'Solvent, blowing agent, gasoline component',
        'molar_mass': 72.15
    },
    'hexane': {
        'iupac_name': 'hexane',
        'common_name': 'hexane',
        'formula': 'C₆H₁₄',
        'smiles': 'CCCCCC',
        'category': 'Alkane',
        'uses': 'Solvent, extraction, adhesives',
        'molar_mass': 86.18
    },
    'cyclohexane': {
        'iupac_name': 'cyclohexane',
        'common_name': 'cyclohexane',
        'formula': 'C₆H₁₂',
        'smiles': 'C1CCCCC1',
        'category': 'Cycloalkane',
        'uses': 'Solvent, nylon production precursor',
        'molar_mass': 84.16
    },
    'phenol': {
        'iupac_name': 'phenol',
        'common_name': 'carbolic acid',
        'formula': 'C₆H₆O',
        'smiles': 'C1=CC=C(C=C1)O',
        'category': 'Phenol',
        'uses': 'Disinfectant, plastic production, antiseptic',
        'molar_mass': 94.11
    },
    'aniline': {
        'iupac_name': 'aniline',
        'common_name': 'aminobenzene',
        'formula': 'C₆H₇N',
        'smiles': 'C1=CC=C(C=C1)N',
        'category': 'Aromatic Amine',
        'uses': 'Dye manufacturing, rubber chemicals',
        'molar_mass': 93.13
    },
    'pyridine': {
        'iupac_name': 'pyridine',
        'common_name': 'pyridine',
        'formula': 'C₅H₅N',
        'smiles': 'C1=CC=NC=C1',
        'category': 'Aromatic Heterocycle',
        'uses': 'Solvent, chemical synthesis, pharmaceuticals',
        'molar_mass': 79.10
    },
    'furan': {
        'iupac_name': 'furan',
        'common_name': 'furan',
        'formula': 'C₄H₄O',
        'smiles': 'C1=COC=C1',
        'category': 'Aromatic Heterocycle',
        'uses': 'Chemical synthesis, polymer production',
        'molar_mass': 68.07
    },
    'dimethyl_sulfoxide': {
        'iupac_name': 'dimethyl sulfoxide',
        'common_name': 'DMSO',
        'formula': 'C₂H₆OS',
        'smiles': 'CS(=O)C',
        'category': 'Organosulfur Compound',
        'uses': 'Solvent, pharmaceutical vehicle, cryoprotectant',
        'molar_mass': 78.13
    },
    'nitromethane': {
        'iupac_name': 'nitromethane',
        'common_name': 'nitromethane',
        'formula': 'CH₃NO₂',
        'smiles': 'C[N+](=O)[O-]',
        'category': 'Nitro Compound',
        'uses': 'Racing fuel, solvent, explosive',
        'molar_mass': 61.04
    },
    'ethylene_oxide': {
        'iupac_name': 'oxirane',
        'common_name': 'ethylene oxide',
        'formula': 'C₂H₄O',
        'smiles': 'C1CO1',
        'category': 'Epoxide',
        'uses': 'Sterilant, chemical intermediate, antifreeze',
        'molar_mass': 44.05
    },
    'propylene_oxide': {
        'iupac_name': '2-methyloxirane',
        'common_name': 'propylene oxide',
        'formula': 'C₃H₆O',
        'smiles': 'CC1CO1',
        'category': 'Epoxide',
        'uses': 'Polyurethane production, fumigant',
        'molar_mass': 58.08
    },
    'trimethylamine': {
        'iupac_name': 'N,N-dimethylmethanamine',
        'common_name': 'trimethylamine',
        'formula': 'C₃H₉N',
        'smiles': 'CN(C)C',
        'category': 'Amine',
        'uses': 'Chemical synthesis, gas treatment',
        'molar_mass': 59.11
    },
    'diethylamine': {
        'iupac_name': 'N-ethylethanamine',
        'common_name': 'diethylamine',
        'formula': 'C₄H₁₁N',
        'smiles': 'CCNCC',
        'category': 'Amine',
        'uses': 'Pharmaceutical synthesis, rubber chemicals',
        'molar_mass': 73.14
    },
    'methylamine': {
        'iupac_name': 'methanamine',
        'common_name': 'methylamine',
        'formula': 'CH₅N',
        'smiles': 'CN',
        'category': 'Amine',
        'uses': 'Chemical synthesis, pharmaceuticals, pesticides',
        'molar_mass': 31.06
    },
    'ethanolamine': {
        'iupac_name': '2-aminoethanol',
        'common_name': 'ethanolamine',
        'formula': 'C₂H₇NO',
        'smiles': 'C(CO)N',
        'category': 'Amine',
        'uses': 'Gas treatment, detergents, cosmetics',
        'molar_mass': 61.08
    },
    'hydrazine': {
        'iupac_name': 'hydrazine',
        'common_name': 'diamine',
        'formula': 'N₂H₄',
        'smiles': 'NN',
        'category': 'Nitrogen Compound',
        'uses': 'Rocket fuel, pharmaceutical synthesis',
        'molar_mass': 32.05
    },
    'hydroxylamine': {
        'iupac_name': 'hydroxylamine',
        'common_name': 'hydroxylamine',
        'formula': 'NH₃O',
        'smiles': 'NO',
        'category': 'Nitrogen Compound',
        'uses': 'Chemical synthesis, reducing agent',
        'molar_mass': 33.03
    },
    'maleic_acid': {
        'iupac_name': '(Z)-but-2-enedioic acid',
        'common_name': 'maleic acid',
        'formula': 'C₄H₄O₄',
        'smiles': 'C(=CC(=O)O)C(=O)O',
        'category': 'Organic Acid',
        'uses': 'Polymer production, food additive',
        'molar_mass': 116.07
    },
    'fumaric_acid': {
        'iupac_name': '(E)-but-2-enedioic acid',
        'common_name': 'fumaric acid',
        'formula': 'C₄H₄O₄',
        'smiles': 'C(=C/C(=O)O)\\C(=O)O',
        'category': 'Organic Acid',
        'uses': 'Food additive, polymer production, pharmaceuticals',
        'molar_mass': 116.07
    },
    'malic_acid': {
        'iupac_name': '2-hydroxybutanedioic acid',
        'common_name': 'malic acid',
        'formula': 'C₄H₆O₅',
        'smiles': 'C(C(C(=O)O)O)C(=O)O',
        'category': 'Organic Acid',
        'uses': 'Food additive, flavor enhancer, cosmetics',
        'molar_mass': 134.09
    },
    'tartaric_acid': {
        'iupac_name': '2,3-dihydroxybutanedioic acid',
        'common_name': 'tartaric acid',
        'formula': 'C₄H₆O₆',
        'smiles': 'C(C(C(=O)O)O)(C(=O)O)O',
        'category': 'Organic Acid',
        'uses': 'Food additive, wine production, pharmaceuticals',
        'molar_mass': 150.09
    },
    
    # Additional High-Value Compounds
    'penicillin': {
        'iupac_name': '(2S,5R,6R)-3,3-dimethyl-7-oxo-6-[(phenylacetyl)amino]-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid',
        'common_name': 'penicillin G',
        'formula': 'C₁₆H₁₈N₂O₄S',
        'smiles': 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
        'category': 'Pharmaceutical',
        'uses': 'Antibiotic, bacterial infection treatment',
        'molar_mass': 334.4
    },
    'morphine': {
        'iupac_name': '(5α,6α)-7,8-didehydro-4,5-epoxy-17-methylmorphinan-3,6-diol',
        'common_name': 'morphine',
        'formula': 'C₁₇H₁₉NO₃',
        'smiles': 'CN1CC[C@]23C4=C5C=CC(=C4[C@H]1CC2)O[C@H]3[C@H](C=C5)O',
        'category': 'Pharmaceutical',
        'uses': 'Pain relief, anesthesia',
        'molar_mass': 285.34
    },
    'testosterone': {
        'iupac_name': '(8R,9S,10R,13S,14S,17S)-17-hydroxy-10,13-dimethyl-1,2,6,7,8,9,10,11,12,13,14,15,16,17-tetradecahydrocyclopenta[a]phenanthren-3-one',
        'common_name': 'testosterone',
        'formula': 'C₁₉H₂₈O₂',
        'smiles': 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C',
        'category': 'Steroid',
        'uses': 'Hormone replacement therapy, muscle building',
        'molar_mass': 288.42
    },
    'adrenaline': {
        'iupac_name': '(R)-4-(1-hydroxy-2-(methylamino)ethyl)benzene-1,2-diol',
        'common_name': 'epinephrine',
        'formula': 'C₉H₁₃NO₃',
        'smiles': 'CNC[C@@H](C1=CC(=C(C=C1)O)O)O',
        'category': 'Pharmaceutical',
        'uses': 'Emergency medicine, allergic reactions',
        'molar_mass': 183.2
    },
    'nicotine': {
        'iupac_name': '(S)-3-(1-methylpyrrolidin-2-yl)pyridine',
        'common_name': 'nicotine',
        'formula': 'C₁₀H₁₄N₂',
        'smiles': 'CN1CCC[C@H]1C2=CN=CC=C2',
        'category': 'Alkaloid',
        'uses': 'Stimulant, insecticide',
        'molar_mass': 162.23
    },
    'adenine': {
        'iupac_name': '9H-purin-6-amine',
        'common_name': 'adenine',
        'formula': 'C₅H₅N₅',
        'smiles': 'C1=NC2=C(N1)C(=NC=N2)N',
        'category': 'Nucleotide Base',
        'uses': 'DNA/RNA component, cellular energy storage',
        'molar_mass': 135.13
    },
    'guanine': {
        'iupac_name': '2-amino-9H-purin-6-ol',
        'common_name': 'guanine',
        'formula': 'C₅H₅N₅O',
        'smiles': 'C1=NC2=C(N1)C(=O)N=C(N2)N',
        'category': 'Nucleotide Base',
        'uses': 'DNA/RNA component',
        'molar_mass': 151.13
    },
    'cytosine': {
        'iupac_name': '4-aminopyrimidin-2-ol',
        'common_name': 'cytosine',
        'formula': 'C₄H₅N₃O',
        'smiles': 'C1=CN=C(N=C1N)O',
        'category': 'Nucleotide Base',
        'uses': 'DNA/RNA component',
        'molar_mass': 111.1
    },
    'thymine': {
        'iupac_name': '5-methylpyrimidine-2,4-diol',
        'common_name': 'thymine',
        'formula': 'C₅H₆N₂O₂',
        'smiles': 'CC1=CN=C(N=C1O)O',
        'category': 'Nucleotide Base',
        'uses': 'DNA component',
        'molar_mass': 126.11
    },
    'capsaicin': {
        'iupac_name': '(E)-N-[(4-hydroxy-3-methoxyphenyl)methyl]-8-methylnon-6-enamide',
        'common_name': 'capsaicin',
        'formula': 'C₁₈H₂₇NO₃',
        'smiles': 'CC(C)CCCC/C=C/C(=O)NCC1=CC(=C(C=C1)O)OC',
        'category': 'Natural Product',
        'uses': 'Spice component, pain relief, pepper spray',
        'molar_mass': 305.41
    },
    'menthol': {
        'iupac_name': '(1R,2S,5R)-2-isopropyl-5-methylcyclohexan-1-ol',
        'common_name': 'menthol',
        'formula': 'C₁₀H₂₀O',
        'smiles': 'C[C@H]1CC[C@@H]([C@H](C1)O)C(C)C',
        'category': 'Alcohol',
        'uses': 'Cooling agent, flavoring, topical analgesic',
        'molar_mass': 156.27
    },
    'vanillin': {
        'iupac_name': '4-hydroxy-3-methoxybenzaldehyde',
        'common_name': 'vanillin',
        'formula': 'C₈H₈O₃',
        'smiles': 'COC1=C(C=CC(=C1)C=O)O',
        'category': 'Aldehyde',
        'uses': 'Flavoring, perfume, food additive',
        'molar_mass': 152.15
    },
    'cinnamaldehyde': {
        'iupac_name': '(E)-3-phenylprop-2-enal',
        'common_name': 'cinnamaldehyde',
        'formula': 'C₉H₈O',
        'smiles': 'C1=CC=C(C=C1)/C=C/C=O',
        'category': 'Aldehyde',
        'uses': 'Cinnamon flavoring, perfume',
        'molar_mass': 132.16
    },
    'resveratrol': {
        'iupac_name': '5-[(E)-2-(4-hydroxyphenyl)ethenyl]benzene-1,3-diol',
        'common_name': 'resveratrol',
        'formula': 'C₁₄H₁₂O₃',
        'smiles': 'C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O',
        'category': 'Polyphenol',
        'uses': 'Antioxidant, dietary supplement',
        'molar_mass': 228.25
    },
    'beta_carotene': {
        'iupac_name': '(1E,3E,5E,7E,9E,11E,13E,15E,17E)-3,7,12,16-tetramethyl-1,18-bis(2,6,6-trimethylcyclohex-1-en-1-yl)octadeca-1,3,5,7,9,11,13,15,17-nonaene',
        'common_name': 'β-carotene',
        'formula': 'C₄₀H₅₆',
        'smiles': 'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(CCCC2(C)C)C)C)C',
        'category': 'Carotenoid',
        'uses': 'Vitamin A precursor, food coloring',
        'molar_mass': 536.87
    },
    'curcumin': {
        'iupac_name': '(1E,6E)-1,7-bis(4-hydroxy-3-methoxyphenyl)hepta-1,6-diene-3,5-dione',
        'common_name': 'curcumin',
        'formula': 'C₂₁H₂₀O₆',
        'smiles': 'COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O',
        'category': 'Polyphenol',
        'uses': 'Spice, anti-inflammatory, antioxidant',
        'molar_mass': 368.38
    },
    'quercetin': {
        'iupac_name': '2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxy-4H-chromen-4-one',
        'common_name': 'quercetin',
        'formula': 'C₁₅H₁₀O₇',
        'smiles': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O',
        'category': 'Flavonoid',
        'uses': 'Antioxidant, dietary supplement',
        'molar_mass': 302.24
    },
    'anthocyanin': {
        'iupac_name': '2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxy-1-benzopyrylium',
        'common_name': 'cyanidin',
        'formula': 'C₁₅H₁₁O₆⁺',
        'smiles': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O',
        'category': 'Flavonoid',
        'uses': 'Natural pigment, antioxidant',
        'molar_mass': 287.24
    },
    'chlorophyll_a': {
        'iupac_name': 'methyl (3S,4S,21S)-14-ethyl-4,8,13,18-tetramethyl-20-oxo-3-(3-oxo-3-{[(2E,7R,11R)-3,7,11,15-tetramethylhexadec-2-en-1-yl]oxy}propyl)-9-vinyl-21,22,23,24-tetradehydroporphyrin-2-carboxylate',
        'common_name': 'chlorophyll a',
        'formula': 'C₅₅H₇₂MgN₄O₅',
        'smiles': 'CC1=C(C2=CC3=NC(=CC4=C(C(=C(N4)C=C5C(=C(C(=N5)C=C(C1=N2)C)C)C=C)C)CCC(=O)OC)C(C3=O)C(=O)OC)C',
        'category': 'Pigment',
        'uses': 'Photosynthesis, natural coloring',
        'molar_mass': 893.49
    },
    'hemoglobin': {
        'iupac_name': '21,22-dihydro-21-(1-hydroxyethyl)-12-methyl-2-vinylporphyrin iron complex',
        'common_name': 'heme',
        'formula': 'C₃₄H₃₂FeN₄O₄',
        'smiles': 'CC1=C(C2=CC3=NC(=CC4=C(C(=C(N4)C=C5C(=C(C(=N5)C=C(C1=N2)C)C)C=C)C)CCC(=O)O)C(C3=O)C(=O)O)C',
        'category': 'Protein Complex',
        'uses': 'Oxygen transport, blood component',
        'molar_mass': 616.49
    }
}

def get_compounds_by_category():
    """Group compounds by category for organized display."""
    categories = {}
    for compound_id, data in COMMON_COMPOUNDS_DATABASE.items():
        category = data['category']
        if category not in categories:
            categories[category] = []
        categories[category].append((compound_id, data))
    return categories

def search_compound(query):
    """Search for compounds by name, formula, or category."""
    query = query.lower()
    results = []
    
    for compound_id, data in COMMON_COMPOUNDS_DATABASE.items():
        if (query in data['iupac_name'].lower() or 
            query in data['common_name'].lower() or
            query in data['formula'].lower() or
            query in data['category'].lower() or
            query in compound_id.lower()):
            results.append((compound_id, data))
    
    return results

def get_total_compounds():
    """Get total number of compounds in database."""
    return len(COMMON_COMPOUNDS_DATABASE)