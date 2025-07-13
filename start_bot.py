#!/usr/bin/env python3
"""
Chemistry Telegram Bot - Universal Starter
Works on any platform with proper error handling and environment detection
"""

import sys
import os
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_dependencies():
    """Check if all required dependencies are available"""
    required_modules = [
        'requests',
        'telegram',
        'PIL',
        'matplotlib',
        'numpy'
    ]
    
    missing = []
    for module in required_modules:
        try:
            __import__(module)
            logger.info(f"✅ {module} - OK")
        except ImportError:
            missing.append(module)
            logger.error(f"❌ {module} - MISSING")
    
    # Check chemistry libraries with better error handling
    try:
        import rdkit
        logger.info("✅ RDKit - OK")
    except ImportError as e:
        logger.warning(f"⚠️ RDKit not available: {e}")
        logger.info("Installing RDKit with: pip install rdkit==2024.3.6")
        # Try to install RDKit automatically
        try:
            import subprocess
            import sys
            subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit==2024.3.6"])
            import rdkit
            logger.info("✅ RDKit installed successfully")
        except Exception as install_error:
            logger.error(f"❌ Failed to install RDKit: {install_error}")
            logger.info("Please install manually: pip install rdkit==2024.3.6")
            missing.append('rdkit')
    
    if missing:
        logger.error(f"Missing dependencies: {', '.join(missing)}")
        logger.info("Install with: pip install -r requirements_deploy.txt")
        return False
    
    return True

def check_environment():
    """Check environment variables and configuration"""
    bot_token = os.getenv('TELEGRAM_BOT_TOKEN')
    if not bot_token:
        logger.error("❌ TELEGRAM_BOT_TOKEN not found in environment")
        logger.info("Set your bot token: export TELEGRAM_BOT_TOKEN=your_token_here")
        return False
    
    logger.info("✅ Bot token configured")
    return True

def main():
    """Main entry point with comprehensive error handling"""
    print("🤖 Chemistry Telegram Bot - Starting...")
    print("=" * 50)
    
    # Check Python version
    if sys.version_info < (3, 8):
        logger.error("❌ Python 3.8+ required")
        sys.exit(1)
    
    logger.info(f"✅ Python {sys.version}")
    
    # Check dependencies
    if not check_dependencies():
        logger.error("❌ Dependency check failed")
        sys.exit(1)
    
    # Check environment
    if not check_environment():
        logger.error("❌ Environment check failed")
        sys.exit(1)
    
    # Import and start bot
    try:
        from telegram_bot_simple import SimpleTelegramBot
        import os
        
        logger.info("🚀 Starting Chemistry Telegram Bot...")
        bot = SimpleTelegramBot(os.getenv('TELEGRAM_BOT_TOKEN'))
        bot.run()
        
    except ImportError as e:
        logger.error(f"❌ Import error: {e}")
        logger.info("Make sure all bot files are in the same directory")
        sys.exit(1)
    except Exception as e:
        logger.error(f"❌ Bot error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n🛑 Bot stopped by user")
        sys.exit(0)
    except Exception as e:
        print(f"❌ Critical error: {e}")
        sys.exit(1)