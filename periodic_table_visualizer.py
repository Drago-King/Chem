"""
Periodic Table Visualizer for Telegram Bot
Creates colorful periodic table images and detailed element information displays.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image, ImageDraw, ImageFont
import io
import numpy as np
from periodic_table_data import PERIODIC_TABLE_DATA, CATEGORY_COLORS

class PeriodicTableVisualizer:
    def __init__(self):
        self.element_data = PERIODIC_TABLE_DATA
        self.category_colors = CATEGORY_COLORS
        
    def create_periodic_table_image(self, width=1200, height=800):
        """Use the high-quality 4K periodic table image provided by user."""
        try:
            # Load the 4K periodic table image
            img_path = 'assets/periodic_table-4k.png'
            img = Image.open(img_path)
            
            # Resize to requested dimensions while maintaining aspect ratio
            img.thumbnail((width, height), Image.Resampling.LANCZOS)
            
            # Create a new image with white background
            new_img = Image.new('RGB', (width, height), 'white')
            
            # Center the periodic table image
            x_offset = (width - img.width) // 2
            y_offset = (height - img.height) // 2
            new_img.paste(img, (x_offset, y_offset))
            
            return new_img
            
        except Exception as e:
            # Fallback to generated table if image fails to load
            return self._create_fallback_periodic_table(width, height)
    
    def _create_fallback_periodic_table(self, width=1200, height=800):
        """Fallback periodic table generator."""
        fig, ax = plt.subplots(figsize=(width/100, height/100))
        ax.set_xlim(0, 18)
        ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        
        # Remove axes
        ax.axis('off')
        
        # Define element positions (group, period)
        element_positions = {
            1: (1, 1), 2: (18, 1),  # Period 1
            3: (1, 2), 4: (2, 2), 5: (13, 2), 6: (14, 2), 7: (15, 2), 8: (16, 2), 9: (17, 2), 10: (18, 2),  # Period 2
            11: (1, 3), 12: (2, 3), 13: (13, 3), 14: (14, 3), 15: (15, 3), 16: (16, 3), 17: (17, 3), 18: (18, 3),  # Period 3
            19: (1, 4), 20: (2, 4), 21: (3, 4), 22: (4, 4), 23: (5, 4), 24: (6, 4), 25: (7, 4), 26: (8, 4), 27: (9, 4), 28: (10, 4), 29: (11, 4), 30: (12, 4), 31: (13, 4), 32: (14, 4), 33: (15, 4), 34: (16, 4), 35: (17, 4), 36: (18, 4),  # Period 4
            37: (1, 5), 38: (2, 5), 39: (3, 5), 40: (4, 5), 41: (5, 5), 42: (6, 5), 43: (7, 5), 44: (8, 5), 45: (9, 5), 46: (10, 5), 47: (11, 5), 48: (12, 5), 49: (13, 5), 50: (14, 5), 51: (15, 5), 52: (16, 5), 53: (17, 5), 54: (18, 5),  # Period 5
            55: (1, 6), 56: (2, 6), 72: (4, 6), 73: (5, 6), 74: (6, 6), 75: (7, 6), 76: (8, 6), 77: (9, 6), 78: (10, 6), 79: (11, 6), 80: (12, 6), 81: (13, 6), 82: (14, 6), 83: (15, 6), 84: (16, 6), 85: (17, 6), 86: (18, 6),  # Period 6
            87: (1, 7), 88: (2, 7), 104: (4, 7), 105: (5, 7), 106: (6, 7), 107: (7, 7), 108: (8, 7), 109: (9, 7), 110: (10, 7), 111: (11, 7), 112: (12, 7), 113: (13, 7), 114: (14, 7), 115: (15, 7), 116: (16, 7), 117: (17, 7), 118: (18, 7),  # Period 7
            # Lanthanides
            57: (3, 8.5), 58: (4, 8.5), 59: (5, 8.5), 60: (6, 8.5), 61: (7, 8.5), 62: (8, 8.5), 63: (9, 8.5), 64: (10, 8.5), 65: (11, 8.5), 66: (12, 8.5), 67: (13, 8.5), 68: (14, 8.5), 69: (15, 8.5), 70: (16, 8.5), 71: (17, 8.5),
            # Actinides
            89: (3, 9.5), 90: (4, 9.5), 91: (5, 9.5), 92: (6, 9.5), 93: (7, 9.5), 94: (8, 9.5), 95: (9, 9.5), 96: (10, 9.5), 97: (11, 9.5), 98: (12, 9.5), 99: (13, 9.5), 100: (14, 9.5), 101: (15, 9.5), 102: (16, 9.5), 103: (17, 9.5)
        }
        
        # Draw elements
        for atomic_num, data in self.element_data.items():
            if atomic_num in element_positions:
                x, y = element_positions[atomic_num]
                x -= 0.4
                y = 10 - y - 0.4  # Flip y-axis
                
                # Get color based on category
                color = data.get('color', '#FFFFFF')
                if 'category' in data:
                    color = self.category_colors.get(data['category'], color)
                
                # Draw element box
                rect = patches.Rectangle((x, y), 0.8, 0.8, linewidth=1, 
                                       edgecolor='black', facecolor=color, alpha=0.8)
                ax.add_patch(rect)
                
                # Add text
                ax.text(x + 0.4, y + 0.6, data['symbol'], ha='center', va='center', 
                       fontsize=8, fontweight='bold')
                ax.text(x + 0.4, y + 0.4, str(atomic_num), ha='center', va='center', 
                       fontsize=6)
                ax.text(x + 0.4, y + 0.2, data['name'][:3], ha='center', va='center', 
                       fontsize=5)
        
        # Add title
        ax.text(9, 9.5, 'Periodic Table of Elements', ha='center', va='center', 
               fontsize=16, fontweight='bold')
        
        # Add legend
        legend_y = 0.5
        for i, (category, color) in enumerate(self.category_colors.items()):
            if i < 5:  # First column
                legend_x = 0.5
                legend_y_pos = legend_y - (i * 0.15)
            else:  # Second column
                legend_x = 9.5
                legend_y_pos = legend_y - ((i - 5) * 0.15)
            
            rect = patches.Rectangle((legend_x, legend_y_pos), 0.3, 0.1, 
                                   facecolor=color, alpha=0.8, edgecolor='black')
            ax.add_patch(rect)
            ax.text(legend_x + 0.4, legend_y_pos + 0.05, category, 
                   ha='left', va='center', fontsize=7)
        
        plt.tight_layout()
        
        # Convert to PIL Image
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img = Image.open(buf)
        plt.close()
        
        return img
    
    def create_element_info_image(self, atomic_number, width=800, height=1000):
        """Create detailed element information image with orbital diagram."""
        if atomic_number not in self.element_data:
            return None
        
        element = self.element_data[atomic_number]
        
        # Create image
        img = Image.new('RGB', (width, height), 'white')
        draw = ImageDraw.Draw(img)
        
        try:
            # Try to use a better font
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 32)
            header_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 20)
            text_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 16)
            small_font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 14)
        except:
            # Fallback to default font
            title_font = ImageFont.load_default()
            header_font = ImageFont.load_default()
            text_font = ImageFont.load_default()
            small_font = ImageFont.load_default()
        
        y_pos = 20
        
        # Title
        title = f"{element['name']} ({element['symbol']})"
        draw.text((width//2, y_pos), title, font=title_font, fill='black', anchor='mt')
        y_pos += 60
        
        # Atomic number and mass
        draw.text((width//2, y_pos), f"Atomic Number: {atomic_number}", 
                 font=header_font, fill='blue', anchor='mt')
        y_pos += 30
        draw.text((width//2, y_pos), f"Atomic Mass: {element['atomic_mass']} u", 
                 font=header_font, fill='blue', anchor='mt')
        y_pos += 50
        
        # Category and color box
        category = element.get('category', 'Unknown')
        color = self.category_colors.get(category, '#CCCCCC')
        
        # Draw category color box
        box_x = width//2 - 100
        box_y = y_pos
        draw.rectangle([box_x, box_y, box_x + 200, box_y + 30], fill=color, outline='black')
        draw.text((width//2, y_pos + 15), category, font=text_font, fill='black', anchor='mm')
        y_pos += 60
        
        # Electron configuration
        draw.text((50, y_pos), "Electron Configuration:", font=header_font, fill='darkgreen')
        y_pos += 25
        draw.text((50, y_pos), element.get('electron_config', 'Unknown'), font=text_font, fill='black')
        y_pos += 40
        
        # Orbital diagram
        draw.text((50, y_pos), "Orbital Diagram:", font=header_font, fill='darkgreen')
        y_pos += 25
        orbital_diagram = element.get('orbital_diagram', 'Complex orbital pattern')
        # Split long orbital diagrams into multiple lines
        if len(orbital_diagram) > 50:
            lines = [orbital_diagram[i:i+50] for i in range(0, len(orbital_diagram), 50)]
            for line in lines:
                draw.text((50, y_pos), line, font=small_font, fill='black')
                y_pos += 20
        else:
            draw.text((50, y_pos), orbital_diagram, font=text_font, fill='black')
            y_pos += 30
        
        y_pos += 20
        
        # Physical properties
        draw.text((50, y_pos), "Physical Properties:", font=header_font, fill='purple')
        y_pos += 30
        
        properties = [
            f"Group: {element.get('group', 'Unknown')}",
            f"Period: {element.get('period', 'Unknown')}",
            f"Melting Point: {element.get('melting_point', 'Unknown')} °C",
            f"Boiling Point: {element.get('boiling_point', 'Unknown')} °C",
            f"Density: {element.get('density', 'Unknown')} g/cm³",
            f"Electronegativity: {element.get('electronegativity', 'Unknown')}"
        ]
        
        for prop in properties:
            draw.text((50, y_pos), prop, font=text_font, fill='black')
            y_pos += 25
        
        # Oxidation states
        y_pos += 10
        draw.text((50, y_pos), "Common Oxidation States:", font=header_font, fill='red')
        y_pos += 25
        ox_states = element.get('oxidation_states', ['Unknown'])
        ox_text = ', '.join(map(str, ox_states))
        draw.text((50, y_pos), ox_text, font=text_font, fill='black')
        y_pos += 40
        
        # Discovery information
        discovery_year = element.get('discovery_year', 'Unknown')
        if discovery_year:
            draw.text((50, y_pos), f"Discovery Year: {discovery_year}", 
                     font=text_font, fill='darkblue')
            y_pos += 30
        
        # Additional interesting facts
        draw.text((50, y_pos), "Interesting Facts:", font=header_font, fill='orange')
        y_pos += 25
        
        # Add element-specific facts
        facts = self._get_element_facts(atomic_number, element)
        for fact in facts:
            # Wrap long facts
            if len(fact) > 60:
                words = fact.split()
                lines = []
                current_line = ""
                for word in words:
                    if len(current_line + word) < 60:
                        current_line += word + " "
                    else:
                        lines.append(current_line.strip())
                        current_line = word + " "
                if current_line:
                    lines.append(current_line.strip())
                
                for line in lines:
                    draw.text((50, y_pos), f"• {line}", font=small_font, fill='black')
                    y_pos += 18
            else:
                draw.text((50, y_pos), f"• {fact}", font=small_font, fill='black')
                y_pos += 20
        
        return img
    
    def _get_element_facts(self, atomic_number, element):
        """Get interesting facts about specific elements."""
        facts_db = {
            1: ["Lightest element in the universe", "Most abundant element in the universe", "Burns with a pale blue flame"],
            2: ["Second most abundant element in the universe", "Used in balloons because it doesn't burn", "Discovered in the Sun before Earth"],
            6: ["Forms more compounds than any other element", "Diamond and graphite are both pure carbon", "Essential for all known life"],
            8: ["Most abundant element in Earth's crust", "Essential for combustion and respiration", "Discovered independently by multiple scientists"],
            26: ["Most common element in Earth's core", "Essential for hemoglobin in blood", "Magnetic at room temperature"],
            79: ["Doesn't tarnish or corrode", "Used as currency for thousands of years", "Excellent conductor of electricity"],
            118: ["Heaviest known element", "Extremely radioactive with very short half-life", "Only a few atoms have ever been created"]
        }
        
        default_facts = [
            f"Belongs to the {element.get('category', 'unknown')} family",
            f"Located in period {element.get('period', '?')} of the periodic table",
            "Has unique chemical and physical properties"
        ]
        
        return facts_db.get(atomic_number, default_facts)