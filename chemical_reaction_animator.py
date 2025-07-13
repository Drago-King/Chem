"""
Chemical Reaction Animation Generator
Creates animated GIFs showing chemical reactions with molecular transformations.
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle, FancyBboxPatch
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import io
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import tempfile
import os

class ChemicalReactionAnimator:
    def __init__(self):
        self.common_reactions = {
            # Combustion Reactions
            'combustion_methane': {
                'equation': 'CH₄ + 2O₂ → CO₂ + 2H₂O',
                'reactants': ['methane', 'oxygen'],
                'products': ['carbon dioxide', 'water'],
                'smiles_reactants': ['C', 'O=O'],
                'smiles_products': ['O=C=O', 'O'],
                'description': 'Methane combustion - burning natural gas',
                'energy_type': 'exothermic'
            },
            'ethanol_combustion': {
                'equation': 'C₂H₅OH + 3O₂ → 2CO₂ + 3H₂O',
                'reactants': ['ethanol', 'oxygen'],
                'products': ['carbon dioxide', 'water'],
                'smiles_reactants': ['CCO', 'O=O'],
                'smiles_products': ['O=C=O', 'O'],
                'description': 'Ethanol combustion - alcohol burning',
                'energy_type': 'exothermic'
            },
            'propane_combustion': {
                'equation': 'C₃H₈ + 5O₂ → 3CO₂ + 4H₂O',
                'reactants': ['propane', 'oxygen'],
                'products': ['carbon dioxide', 'water'],
                'smiles_reactants': ['CCC', 'O=O'],
                'smiles_products': ['O=C=O', 'O'],
                'description': 'Propane combustion - gas grilling fuel',
                'energy_type': 'exothermic'
            },
            'glucose_combustion': {
                'equation': 'C₆H₁₂O₆ + 6O₂ → 6CO₂ + 6H₂O',
                'reactants': ['glucose', 'oxygen'],
                'products': ['carbon dioxide', 'water'],
                'smiles_reactants': ['C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', 'O=O'],
                'smiles_products': ['O=C=O', 'O'],
                'description': 'Glucose combustion - cellular respiration',
                'energy_type': 'exothermic'
            },
            'acetylene_combustion': {
                'equation': '2C₂H₂ + 5O₂ → 4CO₂ + 2H₂O',
                'reactants': ['acetylene', 'oxygen'],
                'products': ['carbon dioxide', 'water'],
                'smiles_reactants': ['C#C', 'O=O'],
                'smiles_products': ['O=C=O', 'O'],
                'description': 'Acetylene combustion - welding torch fuel',
                'energy_type': 'exothermic'
            },
            
            # Synthesis Reactions
            'synthesis_water': {
                'equation': '2H₂ + O₂ → 2H₂O',
                'reactants': ['hydrogen', 'oxygen'],
                'products': ['water'],
                'smiles_reactants': ['[H][H]', 'O=O'],
                'smiles_products': ['O'],
                'description': 'Water synthesis from elements',
                'energy_type': 'exothermic'
            },
            'synthesis_ammonia': {
                'equation': 'N₂ + 3H₂ → 2NH₃',
                'reactants': ['nitrogen', 'hydrogen'],
                'products': ['ammonia'],
                'smiles_reactants': ['N#N', '[H][H]'],
                'smiles_products': ['N'],
                'description': 'Ammonia synthesis - Haber process',
                'energy_type': 'exothermic'
            },
            'synthesis_methanol': {
                'equation': 'CO + 2H₂ → CH₃OH',
                'reactants': ['carbon monoxide', 'hydrogen'],
                'products': ['methanol'],
                'smiles_reactants': ['[C-]#[O+]', '[H][H]'],
                'smiles_products': ['CO'],
                'description': 'Methanol synthesis - industrial alcohol production',
                'energy_type': 'exothermic'
            },
            'synthesis_hydrogen_chloride': {
                'equation': 'H₂ + Cl₂ → 2HCl',
                'reactants': ['hydrogen', 'chlorine'],
                'products': ['hydrogen chloride'],
                'smiles_reactants': ['[H][H]', 'ClCl'],
                'smiles_products': ['Cl'],
                'description': 'Hydrogen chloride synthesis',
                'energy_type': 'exothermic'
            },
            'synthesis_sodium_chloride': {
                'equation': '2Na + Cl₂ → 2NaCl',
                'reactants': ['sodium', 'chlorine'],
                'products': ['sodium chloride'],
                'smiles_reactants': ['[Na]', 'ClCl'],
                'smiles_products': ['[Na+].[Cl-]'],
                'description': 'Salt formation from elements',
                'energy_type': 'exothermic'
            },
            
            # Decomposition Reactions
            'decomposition_water': {
                'equation': '2H₂O → 2H₂ + O₂',
                'reactants': ['water'],
                'products': ['hydrogen', 'oxygen'],
                'smiles_reactants': ['O'],
                'smiles_products': ['[H][H]', 'O=O'],
                'description': 'Water electrolysis - splitting water',
                'energy_type': 'endothermic'
            },
            'decomposition_hydrogen_peroxide': {
                'equation': '2H₂O₂ → 2H₂O + O₂',
                'reactants': ['hydrogen peroxide'],
                'products': ['water', 'oxygen'],
                'smiles_reactants': ['OO'],
                'smiles_products': ['O', 'O=O'],
                'description': 'Hydrogen peroxide decomposition',
                'energy_type': 'exothermic'
            },
            'decomposition_calcium_carbonate': {
                'equation': 'CaCO₃ → CaO + CO₂',
                'reactants': ['calcium carbonate'],
                'products': ['calcium oxide', 'carbon dioxide'],
                'smiles_reactants': ['C(=O)([O-])[O-].[Ca+2]'],
                'smiles_products': ['[Ca+2].[O-2]', 'C(=O)=O'],
                'description': 'Limestone decomposition - cement production',
                'energy_type': 'endothermic'
            },
            'decomposition_ammonium_nitrate': {
                'equation': 'NH₄NO₃ → N₂O + 2H₂O',
                'reactants': ['ammonium nitrate'],
                'products': ['nitrous oxide', 'water'],
                'smiles_reactants': ['[NH4+].[O-][N+](=O)[O-]'],
                'smiles_products': ['N#[N+][O-]', 'O'],
                'description': 'Ammonium nitrate thermal decomposition',
                'energy_type': 'exothermic'
            },
            
            # Acid-Base Reactions
            'acid_base': {
                'equation': 'HCl + NaOH → NaCl + H₂O',
                'reactants': ['hydrochloric acid', 'sodium hydroxide'],
                'products': ['sodium chloride', 'water'],
                'smiles_reactants': ['Cl', '[Na+].[OH-]'],
                'smiles_products': ['[Na+].[Cl-]', 'O'],
                'description': 'Acid-base neutralization',
                'energy_type': 'exothermic'
            },
            'acid_base_sulfuric': {
                'equation': 'H₂SO₄ + 2NaOH → Na₂SO₄ + 2H₂O',
                'reactants': ['sulfuric acid', 'sodium hydroxide'],
                'products': ['sodium sulfate', 'water'],
                'smiles_reactants': ['OS(=O)(=O)O', '[Na+].[OH-]'],
                'smiles_products': ['[Na+].[Na+].[O-]S(=O)(=O)[O-]', 'O'],
                'description': 'Sulfuric acid neutralization',
                'energy_type': 'exothermic'
            },
            'acid_base_acetic': {
                'equation': 'CH₃COOH + NaOH → CH₃COONa + H₂O',
                'reactants': ['acetic acid', 'sodium hydroxide'],
                'products': ['sodium acetate', 'water'],
                'smiles_reactants': ['CC(=O)O', '[Na+].[OH-]'],
                'smiles_products': ['CC(=O)[O-].[Na+]', 'O'],
                'description': 'Acetic acid neutralization',
                'energy_type': 'exothermic'
            },
            'acid_base_carbonic': {
                'equation': 'H₂CO₃ + 2NaOH → Na₂CO₃ + 2H₂O',
                'reactants': ['carbonic acid', 'sodium hydroxide'],
                'products': ['sodium carbonate', 'water'],
                'smiles_reactants': ['C(=O)(O)O', '[Na+].[OH-]'],
                'smiles_products': ['C(=O)([O-])[O-].[Na+].[Na+]', 'O'],
                'description': 'Carbonic acid neutralization',
                'energy_type': 'exothermic'
            },
            
            # Redox Reactions
            'oxidation_iron': {
                'equation': '4Fe + 3O₂ → 2Fe₂O₃',
                'reactants': ['iron', 'oxygen'],
                'products': ['iron oxide'],
                'smiles_reactants': ['[Fe]', 'O=O'],
                'smiles_products': ['[O-2].[O-2].[O-2].[Fe+3].[Fe+3]'],
                'description': 'Iron rusting - oxidation process',
                'energy_type': 'exothermic'
            },
            'reduction_copper_oxide': {
                'equation': 'CuO + H₂ → Cu + H₂O',
                'reactants': ['copper oxide', 'hydrogen'],
                'products': ['copper', 'water'],
                'smiles_reactants': ['[Cu+2].[O-2]', '[H][H]'],
                'smiles_products': ['[Cu]', 'O'],
                'description': 'Copper oxide reduction',
                'energy_type': 'exothermic'
            },
            'zinc_acid_reaction': {
                'equation': 'Zn + 2HCl → ZnCl₂ + H₂',
                'reactants': ['zinc', 'hydrochloric acid'],
                'products': ['zinc chloride', 'hydrogen'],
                'smiles_reactants': ['[Zn]', 'Cl'],
                'smiles_products': ['[Zn+2].[Cl-].[Cl-]', '[H][H]'],
                'description': 'Zinc dissolving in acid - hydrogen production',
                'energy_type': 'exothermic'
            },
            'aluminum_oxide_reduction': {
                'equation': '2Al₂O₃ → 4Al + 3O₂',
                'reactants': ['aluminum oxide'],
                'products': ['aluminum', 'oxygen'],
                'smiles_reactants': ['[Al+3].[Al+3].[O-2].[O-2].[O-2]'],
                'smiles_products': ['[Al]', 'O=O'],
                'description': 'Aluminum smelting - electrolytic reduction',
                'energy_type': 'endothermic'
            },
            
            # Additional Industrial Processes
            'steam_reforming': {
                'equation': 'CH₄ + H₂O → CO + 3H₂',
                'reactants': ['methane', 'steam'],
                'products': ['carbon monoxide', 'hydrogen'],
                'smiles_reactants': ['C', 'O'],
                'smiles_products': ['[C-]#[O+]', '[H][H]'],
                'description': 'Steam reforming for hydrogen production',
                'energy_type': 'endothermic'
            },
            'fischer_tropsch': {
                'equation': 'CO + 2H₂ → CH₂ + H₂O',
                'reactants': ['carbon monoxide', 'hydrogen'],
                'products': ['methylene', 'water'],
                'smiles_reactants': ['[C-]#[O+]', '[H][H]'],
                'smiles_products': ['C', 'O'],
                'description': 'Fischer-Tropsch synthesis',
                'energy_type': 'exothermic'
            },
            'cracking_petroleum': {
                'equation': 'C₁₆H₃₄ → C₈H₁₈ + C₈H₁₆',
                'reactants': ['hexadecane'],
                'products': ['octane', 'octene'],
                'smiles_reactants': ['CCCCCCCCCCCCCCCC'],
                'smiles_products': ['CCCCCCCC', 'CCCCCCCC=C'],
                'description': 'Petroleum cracking for gasoline',
                'energy_type': 'endothermic'
            },
            
            # Pharmaceutical Synthesis
            'aspirin_synthesis': {
                'equation': 'C₇H₆O₃ + C₄H₆O₃ → C₉H₈O₄ + C₂H₄O₂',
                'reactants': ['salicylic acid', 'acetic anhydride'],
                'products': ['aspirin', 'acetic acid'],
                'smiles_reactants': ['C1=CC=C(C(=C1)C(=O)O)O', 'CC(=O)OC(=O)C'],
                'smiles_products': ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(=O)O'],
                'description': 'Aspirin pharmaceutical synthesis',
                'energy_type': 'exothermic'
            },
            'paracetamol_synthesis': {
                'equation': 'C₆H₇NO + C₂H₄O₂ → C₈H₉NO₂ + H₂O',
                'reactants': ['para-aminophenol', 'acetic acid'],
                'products': ['paracetamol', 'water'],
                'smiles_reactants': ['C1=CC(=CC=C1N)O', 'CC(=O)O'],
                'smiles_products': ['CC(=O)NC1=CC=C(C=C1)O', 'O'],
                'description': 'Paracetamol synthesis',
                'energy_type': 'exothermic'
            },
            
            # Environmental Chemistry
            'ozone_formation': {
                'equation': '3O₂ → 2O₃',
                'reactants': ['oxygen'],
                'products': ['ozone'],
                'smiles_reactants': ['O=O'],
                'smiles_products': ['[O-][O+]=O'],
                'description': 'Atmospheric ozone formation',
                'energy_type': 'endothermic'
            },
            'ozone_depletion': {
                'equation': 'O₃ + Cl → ClO + O₂',
                'reactants': ['ozone', 'chlorine radical'],
                'products': ['chlorine monoxide', 'oxygen'],
                'smiles_reactants': ['[O-][O+]=O', '[Cl]'],
                'smiles_products': ['[Cl][O]', 'O=O'],
                'description': 'Ozone layer depletion',
                'energy_type': 'exothermic'
            },
            
            # Food Chemistry
            'maillard_reaction': {
                'equation': 'C₆H₁₂O₆ + NH₃ → complex products',
                'reactants': ['glucose', 'amino compound'],
                'products': ['melanoidins'],
                'smiles_reactants': ['C(C1C(C(C(C(O1)O)O)O)O)O', 'N'],
                'smiles_products': ['C1=CC=C(C=C1)C=O'],  # Simplified product
                'description': 'Browning reaction in cooking',
                'energy_type': 'exothermic'
            },
            'caramelization': {
                'equation': 'C₁₂H₂₂O₁₁ → complex carbohydrate polymers',
                'reactants': ['sucrose'],
                'products': ['caramel'],
                'smiles_reactants': ['C(C1C(C(C(C(O1)O)O)O)O)O'],
                'smiles_products': ['C1=CC=C(C=C1)C=O'],  # Simplified
                'description': 'Sugar caramelization',
                'energy_type': 'endothermic'
            },
            
            # Biological Processes
            'photosynthesis': {
                'equation': '6CO₂ + 6H₂O → C₆H₁₂O₆ + 6O₂',
                'reactants': ['carbon dioxide', 'water'],
                'products': ['glucose', 'oxygen'],
                'smiles_reactants': ['O=C=O', 'O'],
                'smiles_products': ['C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', 'O=O'],
                'description': 'Photosynthesis - plants making glucose',
                'energy_type': 'endothermic'
            },
            'fermentation': {
                'equation': 'C₆H₁₂O₆ → 2C₂H₅OH + 2CO₂',
                'reactants': ['glucose'],
                'products': ['ethanol', 'carbon dioxide'],
                'smiles_reactants': ['C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O'],
                'smiles_products': ['CCO', 'C(=O)=O'],
                'description': 'Alcoholic fermentation - yeast process',
                'energy_type': 'exothermic'
            },
            'protein_hydrolysis': {
                'equation': 'Protein + H₂O → Amino acids',
                'reactants': ['protein', 'water'],
                'products': ['amino acids'],
                'smiles_reactants': ['CC(C(=O)NC(C(=O)O)CC1=CC=CC=C1)N', 'O'],
                'smiles_products': ['C[C@@H](C(=O)O)N', 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N'],
                'description': 'Protein digestion - breaking peptide bonds',
                'energy_type': 'exothermic'
            },
            
            # Industrial Processes
            'steel_production': {
                'equation': 'Fe₂O₃ + 3CO → 2Fe + 3CO₂',
                'reactants': ['iron oxide', 'carbon monoxide'],
                'products': ['iron', 'carbon dioxide'],
                'smiles_reactants': ['[O-2].[O-2].[O-2].[Fe+3].[Fe+3]', '[C-]#[O+]'],
                'smiles_products': ['[Fe]', 'C(=O)=O'],
                'description': 'Steel production - blast furnace reaction',
                'energy_type': 'exothermic'
            },
            'cement_formation': {
                'equation': 'CaCO₃ + SiO₂ → CaSiO₃ + CO₂',
                'reactants': ['calcium carbonate', 'silicon dioxide'],
                'products': ['calcium silicate', 'carbon dioxide'],
                'smiles_reactants': ['C(=O)([O-])[O-].[Ca+2]', 'O=[Si]=O'],
                'smiles_products': ['[Ca+2].[O-][Si]([O-])=O', 'C(=O)=O'],
                'description': 'Cement production - limestone and silica reaction',
                'energy_type': 'endothermic'
            },
            'soap_saponification': {
                'equation': 'C₃H₅(OCOR)₃ + 3NaOH → C₃H₅(OH)₃ + 3RCOONa',
                'reactants': ['triglyceride', 'sodium hydroxide'],
                'products': ['glycerol', 'soap'],
                'smiles_reactants': ['CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC', '[Na+].[OH-]'],
                'smiles_products': ['C(C(CO)O)O', 'CCCCCCCCCCCCCCCC(=O)[O-].[Na+]'],
                'description': 'Soap making - fat saponification',
                'energy_type': 'exothermic'
            }
        }
        
    def create_reaction_animation(self, reaction_name, width=800, height=600, frames=30):
        """Create an animated GIF showing a chemical reaction."""
        if reaction_name not in self.common_reactions:
            return None
            
        reaction = self.common_reactions[reaction_name]
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(width/100, height/100))
        ax.set_xlim(-10, 10)
        ax.set_ylim(-6, 6)
        ax.axis('off')
        ax.set_facecolor('white')
        
        # Animation data storage
        animation_frames = []
        
        # Create frames for the animation with error handling
        try:
            for frame in range(frames):
                ax.clear()
                ax.set_xlim(-10, 10)
                ax.set_ylim(-6, 6)
                ax.axis('off')
                ax.set_facecolor('white')
                
                # Calculate animation progress (0 to 1)
                progress = frame / (frames - 1) if frames > 1 else 0
                
                # Draw reaction elements based on progress
                self._draw_reaction_frame(ax, reaction, progress)
                
                # Save frame to memory with error handling
                try:
                    buf = io.BytesIO()
                    plt.savefig(buf, format='png', bbox_inches='tight', dpi=80, facecolor='white')
                    buf.seek(0)
                    frame_img = Image.open(buf).copy()
                    animation_frames.append(frame_img)
                    buf.close()
                except Exception as e:
                    print(f"Error creating frame {frame}: {e}")
                    continue
        except Exception as e:
            print(f"Error in animation creation: {e}")
            return None
        
        plt.close(fig)
        
        # Create universally compatible animated GIF
        if animation_frames:
            # Convert all frames to standard format for maximum compatibility
            compatible_frames = []
            for frame in animation_frames:
                # Convert to RGB first, then to indexed color for GIF
                if frame.mode != 'RGB':
                    frame = frame.convert('RGB')
                # Resize to standard size
                frame = frame.resize((600, 400), Image.Resampling.LANCZOS)
                # Convert to palette mode with limited colors for compatibility
                frame = frame.quantize(colors=128, method=Image.Quantize.MEDIANCUT)
                compatible_frames.append(frame)
            
            # Create GIF with standard settings for universal compatibility
            gif_buffer = io.BytesIO()
            compatible_frames[0].save(
                gif_buffer,
                format='GIF',
                save_all=True,
                append_images=compatible_frames[1:],
                duration=250,  # Standard timing
                loop=0,
                disposal=2,  # Clear previous frame
                optimize=False  # Disable optimization for compatibility
            )
            gif_buffer.seek(0)
            return gif_buffer
        
        return None
    
    def _draw_reaction_frame(self, ax, reaction, progress):
        """Draw a single frame of the reaction animation."""
        # Title
        ax.text(0, 5, reaction['equation'], ha='center', va='center', 
               fontsize=16, fontweight='bold', color='darkblue')
        
        ax.text(0, 4.2, reaction['description'], ha='center', va='center', 
               fontsize=12, color='gray')
        
        # Draw reactants (left side)
        reactant_x = -6
        reactant_y = 0
        
        # Draw products (right side)
        product_x = 6
        product_y = 0
        
        # Animation phases
        if progress < 0.3:
            # Phase 1: Show reactants clearly
            self._draw_molecules(ax, reactant_x, reactant_y, reaction['reactants'], 'reactants', 1.0)
            self._draw_molecules(ax, product_x, product_y, reaction['products'], 'products', 0.1)
            
        elif progress < 0.7:
            # Phase 2: Transition phase - molecules moving and breaking/forming bonds
            transition_progress = (progress - 0.3) / 0.4
            
            # Reactants fade and move toward center
            reactant_alpha = 1.0 - transition_progress * 0.7
            reactant_pos = reactant_x + transition_progress * 3
            self._draw_molecules(ax, reactant_pos, reactant_y, reaction['reactants'], 'reactants', reactant_alpha)
            
            # Products appear and move from center
            product_alpha = transition_progress * 0.8
            product_pos = product_x - (1 - transition_progress) * 3
            self._draw_molecules(ax, product_pos, product_y, reaction['products'], 'products', product_alpha)
            
            # Draw transition state effects
            self._draw_transition_effects(ax, transition_progress, reaction['energy_type'])
            
        else:
            # Phase 3: Show products clearly
            self._draw_molecules(ax, reactant_x, reactant_y, reaction['reactants'], 'reactants', 0.1)
            self._draw_molecules(ax, product_x, product_y, reaction['products'], 'products', 1.0)
        
        # Draw arrow
        arrow_alpha = min(1.0, progress * 2)
        ax.annotate('', xy=(2, 0), xytext=(-2, 0),
                   arrowprops=dict(arrowstyle='->', lw=3, color='red', alpha=arrow_alpha))
        
        # Energy indicator
        energy_color = 'red' if reaction['energy_type'] == 'exothermic' else 'blue'
        energy_symbol = '+ Heat' if reaction['energy_type'] == 'exothermic' else '+ Energy'
        ax.text(0, -4, energy_symbol, ha='center', va='center', 
               fontsize=12, color=energy_color, alpha=min(1.0, progress))
    
    def _draw_molecules(self, ax, x, y, molecules, mol_type, alpha):
        """Draw molecular representations."""
        colors = {
            'reactants': ['lightblue', 'lightgreen', 'lightyellow'],
            'products': ['lightcoral', 'lightpink', 'lightsalmon']
        }
        
        molecule_colors = colors[mol_type]
        
        for i, molecule in enumerate(molecules):
            mol_x = x
            mol_y = y + (i - len(molecules)/2 + 0.5) * 1.5
            
            # Draw molecule as colored circle with name
            color = molecule_colors[i % len(molecule_colors)]
            circle = Circle((mol_x, mol_y), 0.6, color=color, alpha=alpha, ec='black', linewidth=2)
            ax.add_patch(circle)
            
            # Add molecule name
            ax.text(mol_x, mol_y, molecule[:8], ha='center', va='center', 
                   fontsize=10, fontweight='bold', alpha=alpha)
    
    def _draw_transition_effects(self, ax, progress, energy_type):
        """Draw visual effects during the transition phase."""
        # Draw energy burst effect
        if energy_type == 'exothermic':
            # Red/orange explosion effect
            for i in range(8):
                angle = i * np.pi / 4
                length = progress * 2
                end_x = length * np.cos(angle)
                end_y = length * np.sin(angle)
                ax.plot([0, end_x], [0, end_y], 'r-', alpha=progress, linewidth=3)
        else:
            # Blue energy absorption effect
            for i in range(6):
                angle = i * np.pi / 3
                radius = 1 + progress
                center_x = radius * np.cos(angle)
                center_y = radius * np.sin(angle)
                circle = Circle((center_x, center_y), 0.3, color='blue', alpha=progress*0.5)
                ax.add_patch(circle)
    
    def create_reaction_list_image(self, width=800, height=1000):
        """Create an image showing available reactions."""
        img = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(img)
        
        try:
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 24)
            header_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 18)
            text_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 14)
        except:
            title_font = ImageFont.load_default()
            header_font = ImageFont.load_default()
            text_font = ImageFont.load_default()
        
        y_pos = 30
        
        # Title
        draw.text((width//2, y_pos), "Available Chemical Reactions", 
                 font=title_font, fill='darkblue', anchor='mt')
        y_pos += 60
        
        # Instructions
        draw.text((50, y_pos), "Use /animate [reaction_name] to create animated reaction", 
                 font=header_font, fill='darkgreen')
        y_pos += 50
        
        # List reactions
        for i, (reaction_key, reaction_data) in enumerate(self.common_reactions.items()):
            # Reaction name
            draw.text((50, y_pos), f"{i+1}. {reaction_key.replace('_', ' ').title()}", 
                     font=header_font, fill='purple')
            y_pos += 25
            
            # Equation
            draw.text((70, y_pos), f"Equation: {reaction_data['equation']}", 
                     font=text_font, fill='black')
            y_pos += 20
            
            # Description
            draw.text((70, y_pos), f"Description: {reaction_data['description']}", 
                     font=text_font, fill='gray')
            y_pos += 20
            
            # Energy type
            energy_color = 'red' if reaction_data['energy_type'] == 'exothermic' else 'blue'
            draw.text((70, y_pos), f"Energy: {reaction_data['energy_type'].title()}", 
                     font=text_font, fill=energy_color)
            y_pos += 35
            
            # Command example
            draw.text((70, y_pos), f"Command: /animate {reaction_key}", 
                     font=text_font, fill='darkgreen')
            y_pos += 45
        
        return img
    
    def get_reaction_names(self):
        """Get list of available reaction names."""
        return list(self.common_reactions.keys())
    
    def get_reaction_info(self, reaction_name):
        """Get information about a specific reaction."""
        return self.common_reactions.get(reaction_name)