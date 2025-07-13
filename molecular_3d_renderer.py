import io
import base64
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem

class Molecular3DRenderer:
    """
    A class to render 3D molecular structures as images for Telegram bot display.
    """
    
    def __init__(self):
        # Color mapping for different elements
        self.element_colors = {
            'C': '#000000',  # Carbon - black
            'H': '#FFFFFF',  # Hydrogen - white
            'O': '#FF0000',  # Oxygen - red
            'N': '#0000FF',  # Nitrogen - blue
            'S': '#FFFF00',  # Sulfur - yellow
            'P': '#FFA500',  # Phosphorus - orange
            'Cl': '#00FF00', # Chlorine - green
            'Br': '#A52A2A', # Bromine - brown
            'F': '#90EE90',  # Fluorine - light green
            'I': '#800080',  # Iodine - purple
        }
        
        # Size mapping for different elements
        self.element_sizes = {
            'C': 70,   # Carbon
            'H': 40,   # Hydrogen (smaller)
            'O': 66,   # Oxygen
            'N': 65,   # Nitrogen
            'S': 88,   # Sulfur (larger)
            'P': 98,   # Phosphorus
            'Cl': 79,  # Chlorine
            'Br': 94,  # Bromine
            'F': 57,   # Fluorine
            'I': 115,  # Iodine (largest)
        }
    
    def generate_3d_structure_image(self, mol, width=800, height=600):
        """
        Generate a 3D molecular structure image.
        
        Args:
            mol: RDKit molecule object
            width (int): Image width
            height (int): Image height
            
        Returns:
            PIL.Image: The 3D molecular structure image
        """
        try:
            # Add hydrogens for better 3D representation
            mol_3d = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            
            if result != 0:
                # If 3D embedding fails, try multiple times
                for i in range(10):
                    result = AllChem.EmbedMolecule(mol_3d, randomSeed=i)
                    if result == 0:
                        break
            
            if result == 0:
                # Optimize the 3D structure
                AllChem.MMFFOptimizeMolecule(mol_3d)
                
                # Extract coordinates and create 3D plot
                return self._create_3d_plot(mol_3d, width, height)
            else:
                # Fallback to 2D representation in 3D space
                return self._create_2d_fallback(mol, width, height)
                
        except Exception as e:
            print(f"3D structure generation error: {e}")
            return self._create_2d_fallback(mol, width, height)
    
    def _create_3d_plot(self, mol_3d, width, height):
        """
        Create a 3D matplotlib plot of the molecule.
        
        Args:
            mol_3d: RDKit molecule with 3D coordinates
            width (int): Image width
            height (int): Image height
            
        Returns:
            PIL.Image: The 3D plot as an image
        """
        # Set up the plot
        fig = plt.figure(figsize=(width/100, height/100), dpi=100)
        ax = fig.add_subplot(111, projection='3d')
        
        # Extract atom coordinates
        conformer = mol_3d.GetConformer()
        
        # Plot atoms
        for i, atom in enumerate(mol_3d.GetAtoms()):
            pos = conformer.GetAtomPosition(i)
            element = atom.GetSymbol()
            
            color = self.element_colors.get(element, '#808080')  # Default gray
            size = self.element_sizes.get(element, 60)
            
            ax.scatter([pos.x], [pos.y], [pos.z], 
                      c=color, s=size, alpha=0.8, edgecolors='black', linewidth=0.5)
            
            # Add element labels for non-hydrogen atoms
            if element != 'H':
                ax.text(pos.x, pos.y, pos.z, element, fontsize=8, ha='center', va='center')
        
        # Plot bonds
        for bond in mol_3d.GetBonds():
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            
            pos1 = conformer.GetAtomPosition(atom1_idx)
            pos2 = conformer.GetAtomPosition(atom2_idx)
            
            # Draw bond as a line
            ax.plot([pos1.x, pos2.x], [pos1.y, pos2.y], [pos1.z, pos2.z], 
                   'k-', alpha=0.6, linewidth=1.5)
        
        # Set plot properties
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Molecular Structure')
        
        # Remove axes for cleaner look
        ax.grid(False)
        ax.set_facecolor('white')
        
        # Save to image
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        buf.seek(0)
        
        # Convert to PIL Image
        img = Image.open(buf)
        plt.close(fig)  # Clean up
        
        return img
    
    def _create_2d_fallback(self, mol, width, height):
        """
        Create a 2D representation as a fallback when 3D fails.
        
        Args:
            mol: RDKit molecule object
            width (int): Image width
            height (int): Image height
            
        Returns:
            PIL.Image: The 2D structure image
        """
        try:
            from rdkit.Chem import Draw
            
            # Generate 2D coordinates
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)
            
            # Create 2D image
            img = Draw.MolToImage(mol, size=(width, height))
            return img
            
        except Exception as e:
            # Create a simple error image
            img = Image.new('RGB', (width, height), 'white')
            return img
    
    def create_structure_comparison(self, mol, width=1200, height=600):
        """
        Create a side-by-side comparison of 2D and 3D structures.
        
        Args:
            mol: RDKit molecule object
            width (int): Total image width
            height (int): Image height
            
        Returns:
            PIL.Image: Combined 2D and 3D structure image
        """
        try:
            # Generate 2D structure
            from rdkit.Chem import Draw, rdDepictor
            rdDepictor.Compute2DCoords(mol)
            img_2d = Draw.MolToImage(mol, size=(width//2, height))
            
            # Generate 3D structure
            img_3d = self.generate_3d_structure_image(mol, width//2, height)
            
            # Combine images
            combined = Image.new('RGB', (width, height), 'white')
            combined.paste(img_2d, (0, 0))
            combined.paste(img_3d, (width//2, 0))
            
            return combined
            
        except Exception as e:
            print(f"Structure comparison error: {e}")
            # Return just 2D if comparison fails
            from rdkit.Chem import Draw
            return Draw.MolToImage(mol, size=(width, height))