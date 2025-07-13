import os
import logging
import asyncio
import io
import requests
import json
from chemical_converter import ChemicalConverter
from molecular_3d_renderer import Molecular3DRenderer
from periodic_table_visualizer import PeriodicTableVisualizer
from periodic_table_data import PERIODIC_TABLE_DATA
from chemical_reaction_animator import ChemicalReactionAnimator
from compound_visualizer import CompoundVisualizer
from common_compounds_database import COMMON_COMPOUNDS_DATABASE, get_compounds_by_category, search_compound

# Enable logging
logging.basicConfig(
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO
)
logger = logging.getLogger(__name__)

class SimpleTelegramBot:
    def __init__(self, token):
        self.token = token
        self.base_url = f"https://api.telegram.org/bot{token}"
        self.converter = ChemicalConverter()
        self.renderer_3d = Molecular3DRenderer()
        self.periodic_visualizer = PeriodicTableVisualizer()
        self.reaction_animator = ChemicalReactionAnimator()
        self.compound_visualizer = CompoundVisualizer()
        self.element_data = PERIODIC_TABLE_DATA
        self.compounds_data = COMMON_COMPOUNDS_DATABASE
        self.offset = 0
        
    def send_message(self, chat_id, text, reply_markup=None):
        """Send a text message to a chat."""
        data = {
            'chat_id': chat_id,
            'text': text,
            'parse_mode': 'Markdown'
        }
        if reply_markup:
            data['reply_markup'] = json.dumps(reply_markup)
            
        response = requests.post(f"{self.base_url}/sendMessage", data=data)
        return response.json()
    
    def send_photo(self, chat_id, photo, caption=None, reply_markup=None):
        """Send a photo to a chat."""
        data = {
            'chat_id': chat_id,
            'parse_mode': 'Markdown'
        }
        if caption:
            data['caption'] = caption
        if reply_markup:
            data['reply_markup'] = json.dumps(reply_markup)
            
        files = {'photo': photo}
        response = requests.post(f"{self.base_url}/sendPhoto", data=data, files=files)
        return response.json()
    
    def send_chat_action(self, chat_id, action):
        """Send a chat action (like typing)."""
        data = {
            'chat_id': chat_id,
            'action': action
        }
        response = requests.post(f"{self.base_url}/sendChatAction", data=data)
        return response.json()
    
    def answer_callback_query(self, callback_query_id):
        """Answer a callback query."""
        data = {'callback_query_id': callback_query_id}
        response = requests.post(f"{self.base_url}/answerCallbackQuery", data=data)
        return response.json()
    
    def get_updates(self):
        """Get updates from Telegram."""
        data = {'offset': self.offset, 'timeout': 10}
        response = requests.get(f"{self.base_url}/getUpdates", data=data)
        return response.json()
    
    def handle_start_command(self, chat_id):
        """Handle /start command."""
        welcome_message = """🧪 **Welcome to @IUPAC_name_bot!**
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🎯 **Your Complete Chemistry Assistant**

**Core Features:**
🔬 **Molecular Structures** - Convert IUPAC names to 2D structures
📊 **Periodic Table** - Interactive element exploration with 4K quality
⚗️ **Reaction Animations** - 28+ chemical reaction simulations
📚 **Compound Database** - 104+ compounds across 33 categories
🧬 **Molecular Properties** - Detailed chemical analysis

**Quick Commands:**
• `/table` - View 4K periodic table
• `/ask hydrogen` - Element details with orbital diagrams  
• `/reactions` - Browse 28 reaction animations
• `/animate glucose_combustion` - Specific reaction animation
• `/compounds` - Explore 104+ compound database
• `/compounds pharmaceutical` - Browse by category

**Instant Lookup:**
Just type any chemical name:
• `aspirin` - Get molecular structure
• `caffeine` - See properties & formula
• `glucose` - View 2D structure
• `benzene` - Molecular analysis

**Advanced Features:**
• High-quality molecular visualization
• Comprehensive element data with electron configurations
• Real-time reaction animations
• Searchable compound database
• Educational content for chemistry learning

Ready to explore chemistry? Send any compound name or use the commands above! 🚀

**Animation examples:**
/animate combustion_methane
/animate photosynthesis
/animate acid_base

**Commands:**
/start - Show this help message
/examples - Show quick example buttons
/table - Display periodic table
/ask [element] - Element information
/reactions - Available reaction animations
/animate [reaction] - Create reaction animation
/compounds - Browse compound database
/compounds [category/search] - Find specific compounds

Just send me a chemical name or use the commands!
        """
        self.send_message(chat_id, welcome_message)
    
    def handle_examples_command(self, chat_id):
        """Handle /examples command."""
        keyboard = {
            'inline_keyboard': [
                [
                    {'text': 'Benzene', 'callback_data': 'benzene'},
                    {'text': 'Methane', 'callback_data': 'methane'},
                    {'text': 'Ethanol', 'callback_data': 'ethanol'}
                ],
                [
                    {'text': 'Acetone', 'callback_data': 'acetone'},
                    {'text': 'Glucose', 'callback_data': 'glucose'},
                    {'text': 'Caffeine', 'callback_data': 'caffeine'}
                ],
                [
                    {'text': 'Aspirin', 'callback_data': 'aspirin'},
                    {'text': 'Cholesterol', 'callback_data': 'cholesterol'},
                    {'text': 'Vitamin C', 'callback_data': 'vitamin c'}
                ]
            ]
        }
        
        self.send_message(
            chat_id,
            "🧪 **Quick Examples**\n\nClick any compound below to see its structure:",
            reply_markup=keyboard
        )
    
    def handle_table_command(self, chat_id):
        """Handle /table command - display periodic table."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            "🔬 Generating colorful periodic table..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Generate periodic table image
            table_img = self.periodic_visualizer.create_periodic_table_image()
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            # Convert PIL image to bytes
            bio = io.BytesIO()
            table_img.save(bio, format='PNG')
            bio.seek(0)
            
            caption = """
🧪 **Periodic Table of Elements**

**Color-coded by element categories:**
• Purple: Alkali Metals
• Green: Alkaline Earth Metals  
• Pink: Transition Metals
• Gray: Post-transition Metals
• Orange: Metalloids
• Yellow: Nonmetals
• Light Green: Halogens
• Light Blue: Noble Gases

**How to use:**
Use `/ask [element name]` to get detailed information about any element up to Oganesson (118).

**Examples:**
/ask hydrogen
/ask carbon  
/ask gold
/ask uranium
            """
            
            self.send_photo(chat_id, bio, caption=caption)
            
        except Exception as e:
            logger.error(f"Error generating periodic table: {str(e)}")
            self.send_message(
                chat_id,
                "❌ An error occurred while generating the periodic table. Please try again."
            )
    
    def handle_ask_command(self, chat_id, element_name):
        """Handle /ask command - show detailed element information."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            f"🔬 Looking up detailed information for '{element_name}'..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Find element by name or symbol
            element_found = None
            atomic_number = None
            
            # Search by name (case-insensitive)
            for num, data in self.element_data.items():
                if (data['name'].lower() == element_name.lower() or 
                    data['symbol'].lower() == element_name.lower()):
                    element_found = data
                    atomic_number = num
                    break
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            if element_found:
                # Generate detailed element information image
                element_img = self.periodic_visualizer.create_element_info_image(atomic_number)
                
                if element_img:
                    # Convert PIL image to bytes
                    bio = io.BytesIO()
                    element_img.save(bio, format='PNG')
                    bio.seek(0)
                    
                    caption = f"""
🧪 **{element_found['name']} ({element_found['symbol']})**

**Quick Facts:**
• Atomic Number: {atomic_number}
• Category: {element_found.get('category', 'Unknown')}
• Atomic Mass: {element_found.get('atomic_mass', 'Unknown')} u
• Electron Configuration: {element_found.get('electron_config', 'Unknown')}

**Physical Properties:**
• Group: {element_found.get('group', 'Unknown')}
• Period: {element_found.get('period', 'Unknown')}
• Discovery: {element_found.get('discovery_year', 'Unknown')}

The detailed image shows orbital diagrams, oxidation states, and interesting facts about this element.

Use `/table` to see where this element fits in the periodic table!
                    """
                    
                    self.send_photo(chat_id, bio, caption=caption)
                else:
                    # Fallback text response
                    response_text = f"""
🧪 **{element_found['name']} ({element_found['symbol']})**

**Atomic Information:**
• Atomic Number: {atomic_number}
• Atomic Mass: {element_found.get('atomic_mass', 'Unknown')} u
• Category: {element_found.get('category', 'Unknown')}

**Electron Configuration:**
• {element_found.get('electron_config', 'Unknown')}

**Orbital Diagram:**
• {element_found.get('orbital_diagram', 'Complex orbital pattern')}

**Physical Properties:**
• Group: {element_found.get('group', 'Unknown')}
• Period: {element_found.get('period', 'Unknown')}
• Melting Point: {element_found.get('melting_point', 'Unknown')} °C
• Boiling Point: {element_found.get('boiling_point', 'Unknown')} °C
• Density: {element_found.get('density', 'Unknown')} g/cm³

**Chemical Properties:**
• Electronegativity: {element_found.get('electronegativity', 'Unknown')}
• Oxidation States: {', '.join(map(str, element_found.get('oxidation_states', ['Unknown'])))}

**Discovery:**
• Year: {element_found.get('discovery_year', 'Unknown')}

Use `/table` to see the full periodic table!
                    """
                    
                    self.send_message(chat_id, response_text)
            else:
                # Element not found
                self.send_message(
                    chat_id,
                    f"""
❌ **Element '{element_name}' not found**

**How to use /ask command:**
• `/ask hydrogen` - for Hydrogen
• `/ask au` - for Gold (by symbol)
• `/ask 79` - for Gold (by atomic number)
• `/ask oganesson` - for element 118

**Available elements:** All elements from Hydrogen (1) to Oganesson (118)

Use `/table` to see all available elements in the periodic table!
                    """
                )
                
        except Exception as e:
            logger.error(f"Error processing element request for {element_name}: {str(e)}")
            self.send_message(
                chat_id,
                f"❌ An error occurred while looking up '{element_name}'. Please try again."
            )
    
    def handle_reactions_command(self, chat_id):
        """Handle /reactions command - show available reaction animations."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            "🔬 Loading available chemical reactions..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Generate reactions list image
            reactions_img = self.reaction_animator.create_reaction_list_image()
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            # Convert PIL image to bytes
            bio = io.BytesIO()
            reactions_img.save(bio, format='PNG')
            bio.seek(0)
            
            caption = """
⚗️ **Chemical Reaction Animations**

**Available Reactions:**
1. Combustion Methane
2. Photosynthesis
3. Acid Base Neutralization
4. Water Synthesis
5. Water Decomposition
6. Ethanol Combustion

**How to use:**
Use `/animate [reaction_name]` to create an animated GIF showing the molecular transformation.

**Examples:**
/animate combustion_methane
/animate photosynthesis
/animate acid_base

Each animation shows:
• Reactant molecules on the left
• Product molecules on the right
• Energy changes (exothermic/endothermic)
• Molecular transformation process
            """
            
            self.send_photo(chat_id, bio, caption=caption)
            
        except Exception as e:
            logger.error(f"Error generating reactions list: {str(e)}")
            # Fallback text response
            reactions_text = """
⚗️ **Available Chemical Reaction Animations**

**Reactions you can animate:**

1. **combustion_methane** - CH₄ + 2O₂ → CO₂ + 2H₂O
   Methane burning with oxygen

2. **photosynthesis** - 6CO₂ + 6H₂O → C₆H₁₂O₆ + 6O₂
   Plants making glucose from CO₂ and water

3. **acid_base** - HCl + NaOH → NaCl + H₂O
   Acid and base neutralization

4. **synthesis_water** - 2H₂ + O₂ → 2H₂O
   Making water from hydrogen and oxygen

5. **decomposition_water** - 2H₂O → 2H₂ + O₂
   Splitting water by electrolysis

6. **ethanol_combustion** - C₂H₅OH + 3O₂ → 2CO₂ + 3H₂O
   Alcohol burning in oxygen

**Usage:** `/animate [reaction_name]`
**Example:** `/animate combustion_methane`
            """
            self.send_message(chat_id, reactions_text)
    
    def handle_animate_command(self, chat_id, reaction_name):
        """Handle /animate command - create reaction animation."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            f"🔬 Creating animation for '{reaction_name}' reaction..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Check if reaction exists
            reaction_info = self.reaction_animator.get_reaction_info(reaction_name)
            
            if not reaction_info:
                # Delete processing message
                if processing_msg_id:
                    requests.post(f"{self.base_url}/deleteMessage", data={
                        'chat_id': chat_id,
                        'message_id': processing_msg_id
                    })
                
                available_reactions = self.reaction_animator.get_reaction_names()
                self.send_message(
                    chat_id,
                    f"""
❌ **Reaction '{reaction_name}' not found**

**Available reactions:**
{', '.join(available_reactions)}

**Example usage:**
/animate combustion_methane
/animate photosynthesis

Use `/reactions` to see all available reactions with descriptions.
                    """
                )
                return
            
            # Generate animation with progress feedback
            progress_msg = self.send_message(chat_id, "🎬 Generating animation frames... This may take a moment.")
            progress_msg_id = progress_msg.get('result', {}).get('message_id')
            
            try:
                animation_gif = self.reaction_animator.create_reaction_animation(reaction_name)
            except Exception as anim_error:
                print(f"Animation creation error: {anim_error}")
                animation_gif = None
                
            # Delete progress message
            if progress_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': progress_msg_id
                })
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            if animation_gif:
                # Send the animated GIF with enhanced error handling
                caption = f"""
⚗️ **{reaction_name.replace('_', ' ').title()} Animation**

**Reaction:** {reaction_info['equation']}
**Type:** {reaction_info['energy_type'].title()}
**Description:** {reaction_info['description']}

**Animation shows:**
• Reactants (left) → Products (right)
• Energy changes during reaction
• Molecular transformation process

Use `/reactions` to see more available animations!
                """
                
                try:
                    # Ensure GIF buffer is properly positioned
                    animation_gif.seek(0)
                    
                    # Send as animation (GIF) with proper file naming
                    data = {
                        'chat_id': chat_id,
                        'caption': caption,
                        'parse_mode': 'Markdown'
                    }
                    files = {'animation': (f'{reaction_name}.gif', animation_gif.getvalue(), 'image/gif')}
                    response = requests.post(f"{self.base_url}/sendAnimation", data=data, files=files)
                    
                    if not response.json().get('ok'):
                        # Fallback: send as document with proper MIME type
                        animation_gif.seek(0)
                        files = {'document': (f'{reaction_name}.gif', animation_gif.getvalue(), 'image/gif')}
                        doc_response = requests.post(f"{self.base_url}/sendDocument", data=data, files=files)
                        
                        # If document also fails, try static image as last resort
                        if not doc_response.json().get('ok'):
                            static_image = self._create_static_reaction_image(reaction_info)
                            if static_image:
                                bio = io.BytesIO()
                                static_image.save(bio, format='PNG')
                                bio.seek(0)
                                
                                self.send_photo(chat_id, bio, caption=f"""
⚗️ **{reaction_name.replace('_', ' ').title()} (Static View)**

**Reaction:** {reaction_info['equation']}
**Type:** {reaction_info['energy_type'].title()}
**Description:** {reaction_info['description']}

⚠️ *Animation unavailable - showing static diagram instead*
Use `/reactions` to see all available reactions.""")
                            else:
                                self.send_message(chat_id, f"""⚠️ Animation creation temporarily unavailable.

**Reaction Information:**
• **Name:** {reaction_name.replace('_', ' ').title()}
• **Equation:** {reaction_info['equation']}
• **Type:** {reaction_info['energy_type'].title()}
• **Description:** {reaction_info['description']}

**Troubleshooting:**
• Try again in a few moments
• Use `/reactions` to see all available animations
• Update your Telegram app to the latest version""")
                except Exception as send_error:
                    print(f"Error sending animation: {send_error}")
                    self.send_message(chat_id, f"⚠️ Animation sending failed. Please try again or contact support.")
                
            else:
                self.send_message(
                    chat_id,
                    f"❌ Could not generate animation for '{reaction_name}'. Please try again."
                )
                
        except Exception as e:
            logger.error(f"Error creating animation for {reaction_name}: {str(e)}")
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            self.send_message(
                chat_id,
                f"❌ An error occurred while creating animation for '{reaction_name}'. Please try again."
            )
    
    def handle_compounds_command(self, chat_id, search_query=None):
        """Handle /compounds command - show compound database."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            "🧪 Loading compound database..." if not search_query else f"🧪 Searching for '{search_query}'..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            if not search_query:
                # Show overview of all categories
                compounds_img = self.compound_visualizer.create_compounds_overview_image()
                
                caption = """
🧪 **Common Compounds Database**

**100+ Essential Chemical Compounds**

**Categories included:**
• Organic Acids (acetic acid, citric acid, etc.)
• Alcohols (methanol, ethanol, glycerol, etc.)
• Pharmaceuticals (aspirin, caffeine, ibuprofen, etc.)
• Inorganic Salts (sodium chloride, calcium carbonate, etc.)
• Vitamins (vitamin C, vitamin A, etc.)
• Polymers (polyethylene, nylon, etc.)
• Sugars (glucose, fructose, sucrose, etc.)
• Industrial Chemicals (sulfuric acid, ammonia, etc.)
• And many more categories!

**Usage:**
/compounds [category] - View specific category
/compounds [search term] - Search compounds

**Examples:**
/compounds pharmaceutical
/compounds alcohol
/compounds aspirin
                """
                
            else:
                # Check if it's a category search or compound search
                categories = get_compounds_by_category()
                category_found = False
                
                for category in categories.keys():
                    if search_query.lower() in category.lower():
                        compounds_img = self.compound_visualizer.create_category_compounds_image(search_query)
                        category_found = True
                        
                        caption = f"""
🧪 **{category} Compounds**

**Found {len(categories[category])} compounds in this category**

Each compound shows:
• Common and IUPAC names
• Chemical formula
• Molecular weight
• Primary uses
• SMILES notation

Use compound names directly to get 2D molecular structures!

**Examples:**
Send "aspirin" to see its molecular structure
Send "caffeine" to see its molecular structure
                        """
                        break
                
                if not category_found:
                    # Search for specific compounds
                    search_results = search_compound(search_query)
                    
                    if search_results:
                        compounds_img = self.compound_visualizer.create_search_results_image(search_query)
                        
                        caption = f"""
🧪 **Search Results for '{search_query}'**

**Found {len(search_results)} matching compound(s)**

Each result shows:
• Common and IUPAC names
• Chemical category
• Formula and molecular weight
• Primary uses
• SMILES notation

Send any compound name to get its 2D molecular structure!
                        """
                    else:
                        compounds_img = None
                        
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            if compounds_img:
                # Convert PIL image to bytes
                bio = io.BytesIO()
                compounds_img.save(bio, format='PNG')
                bio.seek(0)
                
                self.send_photo(chat_id, bio, caption=caption)
                
            else:
                # No results found
                available_categories = list(get_compounds_by_category().keys())
                self.send_message(
                    chat_id,
                    f"""
❌ **No compounds found for '{search_query}'**

**Available categories:**
{', '.join(available_categories[:8])}
{', '.join(available_categories[8:16]) if len(available_categories) > 8 else ''}

**Search examples:**
• /compounds organic acid
• /compounds pharmaceutical
• /compounds aspirin
• /compounds alcohol

**Direct compound examples:**
• aspirin
• caffeine
• glucose
• benzene
                    """
                )
                
        except Exception as e:
            logger.error(f"Error processing compounds command: {str(e)}")
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            self.send_message(
                chat_id,
                "❌ An error occurred while loading the compound database. Please try again."
            )
    
    def _format_compound_text(self, compound_key, compound):
        """Format compound information as decorated text."""
        text = f"🧪 **{compound['common_name']}**\n"
        text += f"📋 **IUPAC Name:** {compound['iupac_name']}\n"
        text += f"🧬 **Formula:** {self._enhance_molecular_formula(compound['formula'])}\n"
        text += f"⚖️ **Molar Mass:** {compound['molar_mass']} g/mol\n"
        text += f"🏷️ **Category:** {compound['category']}\n"
        text += f"💼 **Uses:** {compound['uses']}\n"
        if 'smiles' in compound:
            text += f"🔗 **SMILES:** `{compound['smiles']}`"
        return text
    
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
    
    def _create_static_reaction_image(self, reaction_data):
        """Create a static reaction diagram as fallback for animation issues."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            
            fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
            ax.set_xlim(-8, 8)
            ax.set_ylim(-3, 4)
            ax.axis('off')
            
            # Title
            ax.text(0, 3, reaction_data['equation'], ha='center', va='center',
                   fontsize=18, fontweight='bold', color='darkblue')
            
            ax.text(0, 2.3, reaction_data['description'], ha='center', va='center',
                   fontsize=14, color='gray', style='italic')
            
            # Reactants box
            reactant_box = patches.FancyBboxPatch((-7.5, -1), 3, 2, 
                                                boxstyle="round,pad=0.1", 
                                                facecolor='lightblue', 
                                                edgecolor='blue', linewidth=2)
            ax.add_patch(reactant_box)
            ax.text(-6, 0, 'REACTANTS', ha='center', va='center',
                   fontsize=12, fontweight='bold', color='darkblue')
            
            # Arrow
            ax.arrow(-4, 0, 2, 0, head_width=0.3, head_length=0.3, 
                    fc='red', ec='red', linewidth=3)
            ax.text(-3, 0.5, 'REACTION', ha='center', va='center',
                   fontsize=10, fontweight='bold', color='red')
            
            # Products box
            product_box = patches.FancyBboxPatch((4.5, -1), 3, 2,
                                               boxstyle="round,pad=0.1",
                                               facecolor='lightgreen',
                                               edgecolor='green', linewidth=2)
            ax.add_patch(product_box)
            ax.text(6, 0, 'PRODUCTS', ha='center', va='center',
                   fontsize=12, fontweight='bold', color='darkgreen')
            
            # Energy info
            energy_color = 'red' if reaction_data['energy_type'] == 'exothermic' else 'blue'
            ax.text(0, -2.5, f"Energy: {reaction_data['energy_type'].title()}", 
                   ha='center', va='center', fontsize=14, fontweight='bold', 
                   color=energy_color)
            
            # Save to PIL Image
            buffer = BytesIO()
            plt.savefig(buffer, format='PNG', dpi=200, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.close()
            buffer.seek(0)
            
            from PIL import Image
            return Image.open(buffer)
            
        except Exception as e:
            print(f"Error creating static reaction image: {e}")
            return None
    
    def process_chemical_name(self, chat_id, chemical_name):
        """Process a chemical name and send results."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            f"🔬 Converting '{chemical_name}' to molecular structure..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Convert the chemical name
            result = self.converter.convert_iupac_to_structure(chemical_name)
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            if result["success"]:
                # Prepare the response message
                response_text = f"""
🧪 **{chemical_name.title()}**

**Molecular Properties:**
• Formula: `{result['molecular_formula']}`
• Molecular Weight: `{result['molecular_weight']:.2f} g/mol`
• Number of Atoms: `{result['num_atoms']}`
• Number of Bonds: `{result['num_bonds']}`
• SMILES: `{result['smiles']}`

**Additional Properties:**
• LogP: `{result['additional_info']['LogP']}`
• Polar Surface Area: `{result['additional_info']['Polar Surface Area']}`
• Rotatable Bonds: `{result['additional_info']['Rotatable Bonds']}`
• H-Bond Donors: `{result['additional_info']['H-Bond Donors']}`
• H-Bond Acceptors: `{result['additional_info']['H-Bond Acceptors']}`
                """
                
                # Create inline keyboard for 3D visualization
                keyboard = {
                    'inline_keyboard': [
                        [{'text': '🔬 View 3D Structure', 'callback_data': f'3d_{chemical_name}'}]
                    ]
                }
                
                # Send the molecular structure image if available
                if result["structure_image"]:
                    # Convert PIL image to bytes
                    bio = io.BytesIO()
                    result["structure_image"].save(bio, format='PNG')
                    bio.seek(0)
                    
                    self.send_photo(
                        chat_id,
                        bio,
                        caption=response_text,
                        reply_markup=keyboard
                    )
                else:
                    self.send_message(chat_id, response_text, reply_markup=keyboard)
                    
            else:
                # Handle error case
                error_message = f"""
❌ **Conversion Failed**

Could not convert '{chemical_name}' to a molecular structure.

**Possible reasons:**
• Check the spelling of the IUPAC name
• Try a different compound name
• Use standard IUPAC nomenclature

**Suggestions:**
• Try: benzene, methane, ethanol, acetone
• Use /examples for quick test compounds
                """
                
                self.send_message(chat_id, error_message)
                
        except Exception as e:
            logger.error(f"Error processing {chemical_name}: {str(e)}")
            self.send_message(
                chat_id,
                f"❌ An unexpected error occurred while processing '{chemical_name}'. Please try again."
            )
    
    def show_3d_structure(self, chat_id, compound_name):
        """Show 3D molecular structure for a compound."""
        # Send typing action
        self.send_chat_action(chat_id, 'typing')
        
        # Send processing message
        processing_response = self.send_message(
            chat_id,
            f"🔬 Generating 3D structure for '{compound_name}'..."
        )
        processing_msg_id = processing_response.get('result', {}).get('message_id')
        
        try:
            # Convert the chemical name
            result = self.converter.convert_iupac_to_structure(compound_name)
            
            # Delete processing message
            if processing_msg_id:
                requests.post(f"{self.base_url}/deleteMessage", data={
                    'chat_id': chat_id,
                    'message_id': processing_msg_id
                })
            
            if result["success"]:
                # Generate 3D structure using RDKit molecule
                from rdkit import Chem
                mol = Chem.MolFromSmiles(result["smiles"])
                
                if mol:
                    # Generate 3D structure image
                    img_3d = self.renderer_3d.generate_3d_structure_image(mol)
                    
                    if img_3d:
                        # Convert PIL image to bytes
                        bio = io.BytesIO()
                        img_3d.save(bio, format='PNG')
                        bio.seek(0)
                        
                        caption = f"""
🧪 **3D Structure: {compound_name.title()}**

**3D Molecular Visualization**
• Formula: `{result['molecular_formula']}`
• SMILES: `{result['smiles']}`
• Atoms shown with CPK coloring
• Bonds displayed as connecting lines

**Legend:**
• Black: Carbon (C)
• Red: Oxygen (O) 
• Blue: Nitrogen (N)
• White: Hydrogen (H)
• Yellow: Sulfur (S)
                        """
                        
                        self.send_photo(chat_id, bio, caption=caption)
                    else:
                        self.send_message(
                            chat_id,
                            f"❌ Could not generate 3D structure for '{compound_name}'. The molecule may be too complex for 3D visualization."
                        )
                else:
                    self.send_message(
                        chat_id,
                        f"❌ Could not create 3D structure from SMILES notation for '{compound_name}'."
                    )
            else:
                self.send_message(
                    chat_id,
                    f"❌ Could not find compound '{compound_name}' for 3D visualization."
                )
                
        except Exception as e:
            logger.error(f"Error generating 3D structure for {compound_name}: {str(e)}")
            self.send_message(
                chat_id,
                f"❌ An error occurred while generating 3D structure for '{compound_name}'. Please try again."
            )
    
    def handle_callback_query(self, callback_query):
        """Handle callback queries."""
        query_id = callback_query['id']
        chat_id = callback_query['message']['chat']['id']
        data = callback_query['data']
        
        # Answer the callback query
        self.answer_callback_query(query_id)
        
        if data.startswith('3d_'):
            # Handle 3D structure request
            compound_name = data[3:]  # Remove "3d_" prefix
            self.show_3d_structure(chat_id, compound_name)
        else:
            # Handle example compound
            self.process_chemical_name(chat_id, data)
    
    def handle_message(self, message):
        """Handle incoming messages."""
        chat_id = message['chat']['id']
        
        if 'text' in message:
            text = message['text'].strip()
            
            if text == '/start':
                self.handle_start_command(chat_id)
            elif text == '/examples':
                self.handle_examples_command(chat_id)
            elif text == '/table':
                self.handle_table_command(chat_id)
            elif text == '/reactions':
                self.handle_reactions_command(chat_id)
            elif text == '/compounds':
                self.handle_compounds_command(chat_id)
            elif text.startswith('/ask '):
                element_name = text[5:].strip()  # Remove "/ask " prefix
                if element_name:
                    self.handle_ask_command(chat_id, element_name)
                else:
                    self.send_message(chat_id, "Please specify an element name. Example: `/ask hydrogen`")
            elif text.startswith('/animate '):
                reaction_name = text[9:].strip()  # Remove "/animate " prefix
                if reaction_name:
                    self.handle_animate_command(chat_id, reaction_name)
                else:
                    self.send_message(chat_id, "Please specify a reaction name. Example: `/animate combustion_methane`")
            elif text.startswith('/compounds '):
                search_query = text[11:].strip()  # Remove "/compounds " prefix
                if search_query:
                    self.handle_compounds_command(chat_id, search_query)
                else:
                    self.handle_compounds_command(chat_id)
            elif text.startswith('/'):
                self.send_message(chat_id, "Unknown command. Use /start for help.")
            else:
                # Process as chemical name
                self.process_chemical_name(chat_id, text)
    
    def run(self):
        """Run the bot."""
        print("🤖 Chemical Converter Telegram Bot starting...")
        print("Bot is ready to receive messages!")
        
        while True:
            try:
                updates = self.get_updates()
                
                if updates.get('ok'):
                    for update in updates.get('result', []):
                        self.offset = update['update_id'] + 1
                        
                        if 'message' in update:
                            self.handle_message(update['message'])
                        elif 'callback_query' in update:
                            self.handle_callback_query(update['callback_query'])
                
            except Exception as e:
                logger.error(f"Error in main loop: {e}")
                continue

def main():
    """Start the bot."""
    bot_token = os.getenv('TELEGRAM_BOT_TOKEN')
    
    if not bot_token:
        print("Error: TELEGRAM_BOT_TOKEN environment variable not set!")
        return
    
    bot = SimpleTelegramBot(bot_token)
    bot.run()

if __name__ == '__main__':
    main()