"""
Compound Database Visualizer
Creates visual displays for the comprehensive compound database.
"""

from PIL import Image, ImageDraw, ImageFont
import io
from common_compounds_database import COMMON_COMPOUNDS_DATABASE, get_compounds_by_category, search_compound, get_total_compounds

class CompoundVisualizer:
    def __init__(self):
        self.compounds_db = COMMON_COMPOUNDS_DATABASE
        
    def create_compounds_overview_image(self, width=900, height=1200):
        """Create an overview image of all compound categories."""
        img = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(img)
        
        try:
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 28)
            header_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 20)
            text_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 14)
            small_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 12)
        except:
            title_font = ImageFont.load_default()
            header_font = ImageFont.load_default()
            text_font = ImageFont.load_default()
            small_font = ImageFont.load_default()
        
        y_pos = 30
        
        # Title
        total_compounds = get_total_compounds()
        draw.text((width//2, y_pos), f"Common Compounds Database ({total_compounds} compounds)", 
                 font=title_font, fill='darkblue', anchor='mt')
        y_pos += 60
        
        # Instructions
        draw.text((50, y_pos), "Use /compounds [category] or /compounds [search_term]", 
                 font=header_font, fill='darkgreen')
        y_pos += 40
        
        # Get compounds by category
        categories = get_compounds_by_category()
        
        # Category colors
        category_colors = {
            'Organic Acid': '#FF6B6B',
            'Alcohol': '#4ECDC4',
            'Pharmaceutical': '#45B7D1',
            'Inorganic Salt': '#96CEB4',
            'Vitamin': '#FECA57',
            'Polymer': '#FF9FF3',
            'Sugar': '#54A0FF',
            'Inorganic Acid': '#5F27CD',
            'Inorganic Compound': '#00D2D3',
            'Ketone': '#FF6348',
            'Aromatic Hydrocarbon': '#2ED573',
            'Halogenated Solvent': '#FFA502',
            'Steroid': '#3742FA',
            'Amino Acid': '#FF3838',
            'Inorganic Gas': '#70A1FF',
            'Alkane': '#5352ED',
            'Inorganic Base': '#26de81',
            'Inorganic Oxide': '#FC427B'
        }
        
        # Draw categories
        for category, compounds in categories.items():
            color = category_colors.get(category, '#888888')
            
            # Category header with colored background
            header_rect = [40, y_pos - 5, width - 40, y_pos + 25]
            draw.rectangle(header_rect, fill=color, outline='black')
            draw.text((50, y_pos + 10), f"{category} ({len(compounds)} compounds)", 
                     font=header_font, fill='white', anchor='lm')
            y_pos += 40
            
            # List compounds (first 6 per category for better spacing)
            displayed_compounds = compounds[:6]
            for i, (compound_id, data) in enumerate(displayed_compounds):
                if i % 2 == 0:  # Start new row
                    x_left = 60
                    x_right = width // 2 + 20
                    compound_y = y_pos
                else:  # Second column
                    compound_y = y_pos - 20
                
                x = x_left if i % 2 == 0 else x_right
                
                # Enhanced compound name and formula with better subscript rendering
                name = data['common_name'][:22] + "..." if len(data['common_name']) > 22 else data['common_name']
                draw.text((x, compound_y), f"🧪 {name}", font=text_font, fill='black')
                
                # Enhanced molecular formula with proper subscript rendering
                enhanced_formula = self._enhance_molecular_formula(data['formula'])
                draw.text((x + 10, compound_y + 18), f"📋 {enhanced_formula}", font=text_font, fill='darkred')
                
                # Add molar mass for better information
                draw.text((x + 10, compound_y + 35), f"⚖️ {data['molar_mass']} g/mol", font=small_font, fill='gray')
                
                if i % 2 == 1:  # After every two compounds
                    y_pos += 55
            
            if len(displayed_compounds) % 2 == 1:  # If odd number, adjust y_pos
                y_pos += 55
                
            # Show "and X more..." if there are more compounds
            if len(compounds) > 6:
                draw.text((60, y_pos), f"  ... and {len(compounds) - 6} more compounds", 
                         font=small_font, fill='blue')
                y_pos += 25
            
            y_pos += 15  # Space between categories
        
        # Footer instructions
        y_pos += 20
        draw.text((50, y_pos), "Examples:", font=header_font, fill='purple')
        y_pos += 25
        examples = [
            "/compounds organic acid - Show all organic acids",
            "/compounds pharmaceutical - Show all pharmaceuticals", 
            "/compounds aspirin - Search for aspirin",
            "/compounds alcohol - Show all alcohols"
        ]
        
        for example in examples:
            draw.text((70, y_pos), f"• {example}", font=small_font, fill='darkblue')
            y_pos += 18
        
        return img
    
    def _enhance_molecular_formula(self, formula):
        """Enhance molecular formula with better subscript representation."""
        # Convert subscript numbers to more visible format
        subscript_map = {
            '₁': '₁', '₂': '₂', '₃': '₃', '₄': '₄', '₅': '₅', 
            '₆': '₆', '₇': '₇', '₈': '₈', '₉': '₉', '₀': '₀',
            # Also handle regular numbers that should be subscript
            '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅',
            '6': '₆', '7': '₇', '8': '₈', '9': '₉', '0': '₀'
        }
        
        enhanced = ""
        i = 0
        while i < len(formula):
            char = formula[i]
            # Check if current character is a letter followed by numbers
            if char.isalpha():
                enhanced += char
                i += 1
                # Convert following numbers to subscript
                while i < len(formula) and (formula[i].isdigit() or formula[i] in subscript_map):
                    if formula[i].isdigit():
                        enhanced += subscript_map.get(formula[i], formula[i])
                    else:
                        enhanced += formula[i]
                    i += 1
            else:
                enhanced += char
                i += 1
        
        return enhanced
    
    def create_category_compounds_image(self, category, width=800, height=1000):
        """Create detailed view of compounds in a specific category."""
        compounds_by_cat = get_compounds_by_category()
        
        # Find matching category (case-insensitive)
        matching_category = None
        for cat_name in compounds_by_cat.keys():
            if category.lower() in cat_name.lower():
                matching_category = cat_name
                break
        
        if not matching_category:
            return None
            
        compounds = compounds_by_cat[matching_category]
        
        img = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(img)
        
        try:
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 24)
            header_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 18)
            text_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 14)
            small_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 12)
        except:
            title_font = ImageFont.load_default()
            header_font = ImageFont.load_default()
            text_font = ImageFont.load_default()
            small_font = ImageFont.load_default()
        
        y_pos = 30
        
        # Title
        draw.text((width//2, y_pos), f"{matching_category} Compounds", 
                 font=title_font, fill='darkblue', anchor='mt')
        y_pos += 50
        
        draw.text((width//2, y_pos), f"Total: {len(compounds)} compounds", 
                 font=header_font, fill='gray', anchor='mt')
        y_pos += 50
        
        # List compounds with details
        for i, (compound_id, data) in enumerate(compounds):
            if y_pos > height - 150:  # Stop if running out of space
                draw.text((50, y_pos), f"... and {len(compounds) - i} more compounds", 
                         font=text_font, fill='blue')
                break
                
            # Compound box
            box_height = 90
            box_rect = [30, y_pos - 10, width - 30, y_pos + box_height]
            draw.rectangle(box_rect, outline='gray', width=2)
            
            # Compound name
            draw.text((50, y_pos), data['common_name'], font=header_font, fill='darkblue')
            y_pos += 25
            
            # IUPAC name (if different)
            if data['iupac_name'].lower() != data['common_name'].lower():
                iupac_text = f"IUPAC: {data['iupac_name']}"
                if len(iupac_text) > 70:
                    iupac_text = iupac_text[:67] + "..."
                draw.text((50, y_pos), iupac_text, font=small_font, fill='gray')
                y_pos += 18
            
            # Formula and molecular weight
            draw.text((50, y_pos), f"Formula: {data['formula']}", font=text_font, fill='black')
            draw.text((300, y_pos), f"MW: {data.get('molar_mass', 'N/A')} g/mol", font=text_font, fill='black')
            y_pos += 20
            
            # Uses
            uses_text = f"Uses: {data.get('uses', 'Various applications')}"
            if len(uses_text) > 80:
                uses_text = uses_text[:77] + "..."
            draw.text((50, y_pos), uses_text, font=small_font, fill='darkgreen')
            y_pos += 25
        
        return img
    
    def create_search_results_image(self, search_term, width=800, height=1000):
        """Create image showing search results."""
        results = search_compound(search_term)
        
        if not results:
            return None
        
        img = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(img)
        
        try:
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 24)
            header_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 18)
            text_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 14)
            small_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 12)
        except:
            title_font = ImageFont.load_default()
            header_font = ImageFont.load_default()
            text_font = ImageFont.load_default()
            small_font = ImageFont.load_default()
        
        y_pos = 30
        
        # Title
        draw.text((width//2, y_pos), f"Search Results for '{search_term}'", 
                 font=title_font, fill='darkblue', anchor='mt')
        y_pos += 50
        
        draw.text((width//2, y_pos), f"Found {len(results)} compound(s)", 
                 font=header_font, fill='gray', anchor='mt')
        y_pos += 50
        
        # Display results
        for i, (compound_id, data) in enumerate(results):
            if y_pos > height - 150:
                draw.text((50, y_pos), f"... and {len(results) - i} more results", 
                         font=text_font, fill='blue')
                break
                
            # Result box
            box_height = 110
            box_rect = [30, y_pos - 10, width - 30, y_pos + box_height]
            draw.rectangle(box_rect, outline='blue', width=2)
            
            # Compound name and category
            draw.text((50, y_pos), data['common_name'], font=header_font, fill='darkblue')
            draw.text((width - 200, y_pos), f"[{data['category']}]", font=text_font, fill='purple')
            y_pos += 25
            
            # IUPAC name
            if data['iupac_name'].lower() != data['common_name'].lower():
                iupac_text = f"IUPAC: {data['iupac_name']}"
                if len(iupac_text) > 70:
                    iupac_text = iupac_text[:67] + "..."
                draw.text((50, y_pos), iupac_text, font=small_font, fill='gray')
                y_pos += 18
            
            # Formula and molecular weight
            draw.text((50, y_pos), f"Formula: {data['formula']}", font=text_font, fill='black')
            draw.text((300, y_pos), f"MW: {data.get('molar_mass', 'N/A')} g/mol", font=text_font, fill='black')
            y_pos += 20
            
            # Uses
            uses_text = f"Uses: {data.get('uses', 'Various applications')}"
            if len(uses_text) > 80:
                uses_text = uses_text[:77] + "..."
            draw.text((50, y_pos), uses_text, font=small_font, fill='darkgreen')
            y_pos += 25
            
            # SMILES notation
            smiles_text = f"SMILES: {data.get('smiles', 'N/A')}"
            if len(smiles_text) > 70:
                smiles_text = smiles_text[:67] + "..."
            draw.text((50, y_pos), smiles_text, font=small_font, fill='orange')
            y_pos += 30
        
        return img