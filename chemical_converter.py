from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
import base64
import requests
import json
import time

class ChemicalConverter:
    """
    A class to convert IUPAC chemical names to molecular structures and properties.
    """
    
    def __init__(self):
        """Initialize the chemical converter."""
        pass
    
    def convert_iupac_to_structure(self, iupac_name):
        """
        Convert IUPAC name to molecular structure and properties.
        
        Args:
            iupac_name (str): The IUPAC name of the chemical compound
            
        Returns:
            dict: A dictionary containing structure data and properties
        """
        try:
            # Clean the input
            iupac_name = iupac_name.strip()
            
            if not iupac_name:
                return {
                    "success": False,
                    "error": "Please enter a valid IUPAC name"
                }
            
            # Try to create molecule from IUPAC name
            mol = None
            
            # First, check if it's in our common compounds dictionary
            mol = self._get_molecule_from_common_names(iupac_name)
            
            # If not found, try direct conversion from name
            if mol is None:
                try:
                    mol = Chem.MolFromName(iupac_name)
                except:
                    mol = None
            
            # If that fails, try some common name variations
            if mol is None:
                # Try with different capitalizations and common variations
                variations = [
                    iupac_name.lower(),
                    iupac_name.upper(),
                    iupac_name.capitalize(),
                    iupac_name.replace('-', ''),
                    iupac_name.replace(' ', ''),
                ]
                
                for variation in variations:
                    try:
                        mol = Chem.MolFromName(variation)
                        if mol is not None:
                            break
                    except:
                        continue
                    
                    # Also try our common names dictionary for variations
                    mol = self._get_molecule_from_common_names(variation)
                    if mol is not None:
                        break
            
            # If still not found, try multiple external databases
            if mol is None:
                # Try PubChem API for complex IUPAC names
                mol = self._get_molecule_from_pubchem(iupac_name)
            
            # Try ChEBI database if PubChem fails
            if mol is None:
                mol = self._get_molecule_from_chebi(iupac_name)
            
            # Try NCI/CADD database
            if mol is None:
                mol = self._get_molecule_from_nci(iupac_name)
            
            # Try alternative PubChem search methods
            if mol is None:
                mol = self._get_molecule_pubchem_alternative(iupac_name)
            
            if mol is None:
                return {
                    "success": False,
                    "error": f"Could not convert '{iupac_name}' to a valid molecular structure. "
                            "Please check the IUPAC name spelling or try a different compound."
                }
            
            # Generate 2D coordinates for drawing
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)
            
            # Generate molecular structure image
            structure_image = self._generate_structure_image(mol)
            
            # Generate 3D structure data
            structure_3d_data = self._generate_3d_structure(mol)
            
            # Calculate molecular properties
            properties = self._calculate_properties(mol)
            
            # Get SMILES notation
            smiles = Chem.MolToSmiles(mol)
            
            return {
                "success": True,
                "structure_image": structure_image,
                "structure_3d": structure_3d_data,
                "molecular_formula": properties["formula"],
                "molecular_weight": properties["molecular_weight"],
                "num_atoms": properties["num_atoms"],
                "num_bonds": properties["num_bonds"],
                "smiles": smiles,
                "additional_info": {
                    "LogP": f"{properties['logp']:.2f}",
                    "Polar Surface Area": f"{properties['tpsa']:.2f} Ų",
                    "Rotatable Bonds": properties["rotatable_bonds"],
                    "H-Bond Donors": properties["hbd"],
                    "H-Bond Acceptors": properties["hba"]
                }
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": f"An error occurred while processing '{iupac_name}': {str(e)}"
            }
    
    def _generate_structure_image(self, mol, width=400, height=300):
        """
        Generate a 2D structure image of the molecule.
        
        Args:
            mol: RDKit molecule object
            width (int): Image width
            height (int): Image height
            
        Returns:
            PIL.Image: The molecular structure image
        """
        try:
            # Create a drawer
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            
            # Set drawing options for better visualization
            opts = drawer.drawOptions()
            opts.addStereoAnnotation = True
            opts.addAtomIndices = False
            opts.bondLineWidth = 2
            
            # Draw the molecule
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            
            # Get the image data
            img_data = drawer.GetDrawingText()
            
            # Convert to PIL Image
            img = Image.open(io.BytesIO(img_data))
            
            return img
            
        except Exception as e:
            # Fallback to basic RDKit drawing if Cairo fails
            try:
                img = Draw.MolToImage(mol, size=(width, height))
                return img
            except Exception as fallback_e:
                # Return None if both methods fail
                return None
    
    def _calculate_properties(self, mol):
        """
        Calculate various molecular properties.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            dict: Dictionary containing molecular properties
        """
        try:
            # Basic properties
            formula = rdMolDescriptors.CalcMolFormula(mol)
            molecular_weight = Descriptors.MolWt(mol)
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            
            # Advanced properties
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            return {
                "formula": formula,
                "molecular_weight": molecular_weight,
                "num_atoms": num_atoms,
                "num_bonds": num_bonds,
                "logp": logp,
                "tpsa": tpsa,
                "rotatable_bonds": rotatable_bonds,
                "hbd": hbd,
                "hba": hba
            }
            
        except Exception as e:
            # Return basic properties if advanced calculation fails
            return {
                "formula": "Unknown",
                "molecular_weight": 0.0,
                "num_atoms": mol.GetNumAtoms() if mol else 0,
                "num_bonds": mol.GetNumBonds() if mol else 0,
                "logp": 0.0,
                "tpsa": 0.0,
                "rotatable_bonds": 0,
                "hbd": 0,
                "hba": 0
            }
    
    def _get_molecule_from_common_names(self, name):
        """
        Get molecule from a dictionary of common chemical names and their SMILES.
        
        Args:
            name (str): Chemical name to look up
            
        Returns:
            mol: RDKit molecule object or None if not found
        """
        # Dictionary of common compounds with their SMILES notation
        common_compounds = {
            # Basic alkanes
            "methane": "C",
            "ethane": "CC",
            "propane": "CCC",
            "butane": "CCCC",
            "pentane": "CCCCC",
            "hexane": "CCCCCC",
            "heptane": "CCCCCCC",
            "octane": "CCCCCCCC",
            
            # Branched alkanes
            "2-methylpropane": "CC(C)C",
            "isobutane": "CC(C)C",
            "2-methylbutane": "CC(C)CC",
            "isopentane": "CC(C)CC",
            "2,2-dimethylpropane": "CC(C)(C)C",
            "neopentane": "CC(C)(C)C",
            
            # Alkenes
            "ethene": "C=C",
            "ethylene": "C=C",
            "propene": "CC=C",
            "propylene": "CC=C",
            "1-butene": "CCC=C",
            "2-butene": "CC=CC",
            
            # Alkynes
            "ethyne": "C#C",
            "acetylene": "C#C",
            "propyne": "CC#C",
            "1-butyne": "CCC#C",
            "2-butyne": "CC#CC",
            
            # Aromatics
            "benzene": "c1ccccc1",
            "toluene": "Cc1ccccc1",
            "methylbenzene": "Cc1ccccc1",
            "xylene": "Cc1ccccc1C",
            "phenol": "Oc1ccccc1",
            "aniline": "Nc1ccccc1",
            "styrene": "C=Cc1ccccc1",
            "naphthalene": "c1ccc2ccccc2c1",
            
            # Alcohols
            "methanol": "CO",
            "ethanol": "CCO",
            "propanol": "CCCO",
            "1-propanol": "CCCO",
            "2-propanol": "CC(C)O",
            "isopropanol": "CC(C)O",
            "butanol": "CCCCO",
            "1-butanol": "CCCCO",
            "2-butanol": "CC(O)CC",
            "2-methyl-2-propanol": "CC(C)(C)O",
            "tert-butanol": "CC(C)(C)O",
            
            # Ethers
            "diethyl ether": "CCOCC",
            "ether": "CCOCC",
            "dimethyl ether": "COC",
            "methyl tert-butyl ether": "COC(C)(C)C",
            "mtbe": "COC(C)(C)C",
            
            # Aldehydes
            "formaldehyde": "C=O",
            "acetaldehyde": "CC=O",
            "ethanal": "CC=O",
            "propanal": "CCC=O",
            "propionaldehyde": "CCC=O",
            "butanal": "CCCC=O",
            "benzaldehyde": "O=Cc1ccccc1",
            
            # Ketones
            "acetone": "CC(=O)C",
            "propanone": "CC(=O)C",
            "butanone": "CCC(=O)C",
            "methyl ethyl ketone": "CCC(=O)C",
            "2-pentanone": "CCCC(=O)C",
            "cyclohexanone": "O=C1CCCCC1",
            
            # Carboxylic acids
            "formic acid": "C(=O)O",
            "methanoic acid": "C(=O)O",
            "acetic acid": "CC(=O)O",
            "ethanoic acid": "CC(=O)O",
            "propionic acid": "CCC(=O)O",
            "propanoic acid": "CCC(=O)O",
            "butyric acid": "CCCC(=O)O",
            "butanoic acid": "CCCC(=O)O",
            "benzoic acid": "O=C(O)c1ccccc1",
            
            # Esters
            "methyl acetate": "CC(=O)OC",
            "ethyl acetate": "CC(=O)OCC",
            "methyl formate": "COC=O",
            "ethyl formate": "CCOC=O",
            "methyl propanoate": "CCC(=O)OC",
            "ethyl propanoate": "CCC(=O)OCC",
            
            # Amines
            "methylamine": "CN",
            "ethylamine": "CCN",
            "propylamine": "CCCN",
            "dimethylamine": "CNC",
            "diethylamine": "CCNCC",
            "trimethylamine": "CN(C)C",
            "triethylamine": "CCN(CC)CC",
            
            # Sugars
            "glucose": "OCC1OC(O)C(O)C(O)C1O",
            "fructose": "OCC1(O)OCC(O)C(O)C1O",
            "sucrose": "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
            
            # Common molecules
            "water": "O",
            "ammonia": "N",
            "hydrogen sulfide": "S",
            "carbon dioxide": "O=C=O",
            "carbon monoxide": "C#O",
            "hydrogen peroxide": "OO",
            "sulfuric acid": "OS(=O)(=O)O",
            "nitric acid": "O[N+](=O)[O-]",
            "hydrochloric acid": "Cl",
            
            # Pharmaceuticals and common compounds
            "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
            "acetylsalicylic acid": "CC(=O)Oc1ccccc1C(=O)O",
            "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "nicotine": "CN1CCCC1c2cccnc2",
            "menthol": "CC(C)C1CCC(C)CC1O",
            "vanillin": "COc1cc(C=O)ccc1O",
            "citric acid": "OC(=O)CC(O)(C(=O)O)CC(=O)O",
            
            # Polymers precursors
            "vinyl chloride": "ClC=C",
            "chloroethene": "ClC=C",
            "acrylonitrile": "C=CC#N",
            "propylene glycol": "CC(O)CO",
            "ethylene glycol": "OCCO",
            
            # Additional complex compounds
            "cyclohexane": "C1CCCCC1",
            "cyclopentane": "C1CCCC1",
            "cyclobutane": "C1CCC1",
            "cyclopropane": "C1CC1",
            "methylcyclohexane": "CC1CCCCC1",
            "ethylbenzene": "CCc1ccccc1",
            "cumene": "CC(C)c1ccccc1",
            "isopropylbenzene": "CC(C)c1ccccc1",
            
            # Heterocycles
            "pyridine": "c1ccncc1",
            "pyrrole": "c1cc[nH]c1",
            "furane": "c1ccoc1",
            "thiophene": "c1ccsc1",
            "imidazole": "c1c[nH]cn1",
            "pyrazole": "c1cn[nH]c1",
            "oxazole": "c1coc[nH]1",
            "thiazole": "c1csc[nH]1",
            
            # Amino acids
            "glycine": "NCC(=O)O",
            "alanine": "CC(N)C(=O)O",
            "valine": "CC(C)C(N)C(=O)O",
            "leucine": "CC(C)CC(N)C(=O)O",
            "isoleucine": "CCC(C)C(N)C(=O)O",
            "serine": "OCC(N)C(=O)O",
            "threonine": "CC(O)C(N)C(=O)O",
            "tyrosine": "Oc1ccc(CC(N)C(=O)O)cc1",
            "phenylalanine": "c1ccc(CC(N)C(=O)O)cc1",
            
            # Nucleotides and bases
            "adenine": "c1nc(N)c2ncnc2n1",
            "guanine": "NC1=Nc2c(ncn2)C(=O)N1",
            "cytosine": "NC1=CC(=O)NC=N1",
            "thymine": "CC1=CNC(=O)NC1=O",
            "uracil": "O=C1C=CNC(=O)N1",
            
            # Steroids
            "cholesterol": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
            "testosterone": "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",
            "estradiol": "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O",
            
            # Natural products
            "thymol": "Cc1ccc(C(C)C)c(O)c1",
            "carvone": "CC(C)C1CC=C(C)C(=O)C1",
            "limonene": "CC(C)C1CC=C(C)CC1",
            "pinene": "CC1=CCC2CC1C2(C)C",
            
            # Vitamins
            "ascorbic acid": "OC1=C(O)C(=O)OC1C(O)CO",
            "vitamin c": "OC1=C(O)C(=O)OC1C(O)CO",
            "thiamine": "OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc2N",
            "riboflavin": "Cc1cc2nc3c(=O)[nH]c(=O)nc3n(CC(O)C(O)C(O)CO)c2cc1C",
            
            # Dyes and pigments
            "indigo": "O=C1/C(=C\\c2[nH]c3ccccc3c2=O)Nc2ccccc12",
            "alizarin": "O=C1c2ccccc2C(=O)c2c(O)ccc(O)c12",
            "methylene blue": "CN(C)c1ccc2nc3ccc(N(C)C)cc3[s+]c2c1",
            
            # Pesticides and herbicides
            "ddt": "ClC(Cl)C(c1ccc(Cl)cc1)c2ccc(Cl)cc2",
            "malathion": "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC",
            "parathion": "CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1",
        }
        
        # Normalize the input name for comparison
        normalized_name = name.lower().strip()
        
        # Direct lookup
        if normalized_name in common_compounds:
            try:
                return Chem.MolFromSmiles(common_compounds[normalized_name])
            except:
                return None
        
        # Try without hyphens and spaces
        normalized_no_punct = normalized_name.replace('-', '').replace(' ', '')
        for key, smiles in common_compounds.items():
            if key.replace('-', '').replace(' ', '') == normalized_no_punct:
                try:
                    return Chem.MolFromSmiles(smiles)
                except:
                    continue
        
        return None
    
    def _get_molecule_from_pubchem(self, name):
        """
        Get molecule from PubChem API using IUPAC name.
        
        Args:
            name (str): Chemical name to look up
            
        Returns:
            mol: RDKit molecule object or None if not found
        """
        try:
            # Clean the name for API request
            name_clean = name.strip()
            
            # PubChem REST API endpoint for name to SMILES conversion
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_clean}/property/IsomericSMILES/JSON"
            
            # Make request with timeout
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and len(properties) > 0:
                        smiles = properties[0].get('IsomericSMILES')
                        if smiles:
                            try:
                                mol = Chem.MolFromSmiles(smiles)
                                return mol
                            except:
                                return None
            
            # If IsomericSMILES fails, try CanonicalSMILES
            url_canonical = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_clean}/property/CanonicalSMILES/JSON"
            response = requests.get(url_canonical, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and len(properties) > 0:
                        smiles = properties[0].get('CanonicalSMILES')
                        if smiles:
                            try:
                                mol = Chem.MolFromSmiles(smiles)
                                return mol
                            except:
                                return None
            
            return None
            
        except Exception as e:
            # If API fails, return None to continue with other methods
            return None
    
    def _generate_3d_structure(self, mol):
        """
        Generate 3D structure data for molecular visualization.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            dict: 3D structure data including SDF and coordinates
        """
        try:
            from rdkit.Chem import AllChem
            
            # Add hydrogens for 3D structure
            mol_3d = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            
            if result == 0:  # Success
                # Optimize the 3D structure
                AllChem.MMFFOptimizeMolecule(mol_3d)
                
                # Generate SDF format for 3D visualization
                sdf_data = Chem.MolToMolBlock(mol_3d)
                
                # Extract atom coordinates
                conformer = mol_3d.GetConformer()
                atoms = []
                
                for i, atom in enumerate(mol_3d.GetAtoms()):
                    pos = conformer.GetAtomPosition(i)
                    atoms.append({
                        'element': atom.GetSymbol(),
                        'x': pos.x,
                        'y': pos.y,
                        'z': pos.z
                    })
                
                # Extract bonds
                bonds = []
                for bond in mol_3d.GetBonds():
                    bonds.append({
                        'atom1': bond.GetBeginAtomIdx(),
                        'atom2': bond.GetEndAtomIdx(),
                        'order': bond.GetBondType().name
                    })
                
                return {
                    'sdf': sdf_data,
                    'atoms': atoms,
                    'bonds': bonds,
                    'has_3d': True
                }
            else:
                return {
                    'sdf': None,
                    'atoms': [],
                    'bonds': [],
                    'has_3d': False,
                    'error': '3D structure generation failed'
                }
                
        except Exception as e:
            return {
                'sdf': None,
                'atoms': [],
                'bonds': [],
                'has_3d': False,
                'error': f'3D generation error: {str(e)}'
            }
    
    def _get_molecule_from_chebi(self, name):
        """
        Get molecule from ChEBI database using name.
        
        Args:
            name (str): Chemical name to look up
            
        Returns:
            mol: RDKit molecule object or None if not found
        """
        try:
            # ChEBI REST API endpoint
            url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{name}"
            
            # Try direct search first
            search_url = f"https://www.ebi.ac.uk/webservices/chebi/2.0/test/getLiteEntity?search={name}&searchCategory=ALL&maxResults=1&stars=3"
            
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200:
                # Parse response for SMILES
                if "smiles" in response.text.lower():
                    # Simple parsing for SMILES notation
                    lines = response.text.split('\n')
                    for line in lines:
                        if 'smiles' in line.lower():
                            smiles = line.split('>')[-1].split('<')[0].strip()
                            if smiles and len(smiles) > 2:
                                try:
                                    mol = Chem.MolFromSmiles(smiles)
                                    if mol:
                                        return mol
                                except:
                                    continue
            
            return None
            
        except Exception as e:
            return None
    
    def _get_molecule_from_nci(self, name):
        """
        Get molecule from NCI/CADD database.
        
        Args:
            name (str): Chemical name to look up
            
        Returns:
            mol: RDKit molecule object or None if not found
        """
        try:
            # NCI CADD Chemical Identifier Resolver
            url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
            
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                smiles = response.text.strip()
                if smiles and not smiles.startswith("Error") and len(smiles) > 1:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        return mol
                    except:
                        return None
            
            return None
            
        except Exception as e:
            return None
    
    def _get_molecule_pubchem_alternative(self, name):
        """
        Alternative PubChem search methods for better coverage.
        
        Args:
            name (str): Chemical name to look up
            
        Returns:
            mol: RDKit molecule object or None if not found
        """
        try:
            # Try synonym search
            synonym_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/synonyms/JSON"
            
            response = requests.get(synonym_url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'InformationList' in data and 'Information' in data['InformationList']:
                    info = data['InformationList']['Information']
                    if info and len(info) > 0:
                        cid = info[0].get('CID')
                        if cid:
                            # Get SMILES using CID
                            smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON"
                            smiles_response = requests.get(smiles_url, timeout=10)
                            
                            if smiles_response.status_code == 200:
                                smiles_data = smiles_response.json()
                                if 'PropertyTable' in smiles_data and 'Properties' in smiles_data['PropertyTable']:
                                    properties = smiles_data['PropertyTable']['Properties']
                                    if properties and len(properties) > 0:
                                        smiles = properties[0].get('IsomericSMILES')
                                        if smiles:
                                            try:
                                                mol = Chem.MolFromSmiles(smiles)
                                                return mol
                                            except:
                                                return None
            
            # Try formula-based search if name might be a formula
            if name.replace(' ', '').replace('-', '').replace('(', '').replace(')', '').isalnum():
                formula_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{name}/property/IsomericSMILES/JSON"
                response = requests.get(formula_url, timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                        properties = data['PropertyTable']['Properties']
                        if properties and len(properties) > 0:
                            smiles = properties[0].get('IsomericSMILES')
                            if smiles:
                                try:
                                    mol = Chem.MolFromSmiles(smiles)
                                    return mol
                                except:
                                    return None
            
            return None
            
        except Exception as e:
            return None
