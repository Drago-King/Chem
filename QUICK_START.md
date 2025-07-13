# 🚀 Chemistry Telegram Bot - Quick Start Guide

## 📱 Start Commands

### Option 1: Enhanced Starter (Recommended)
```bash
python3 start_bot.py
```

### Option 2: Direct Start
```bash
python3 telegram_bot_simple.py
```

### Option 3: Run Script
```bash
python3 run.py
```

## 🔧 Environment Setup

1. **Set your bot token:**
   ```bash
   export TELEGRAM_BOT_TOKEN=your_bot_token_here
   ```

2. **Or create .env file:**
   ```bash
   echo "TELEGRAM_BOT_TOKEN=your_bot_token_here" > .env
   ```

## ✅ Bot Status Check

The enhanced starter will automatically check:
- ✅ Python version (3.8+)
- ✅ All dependencies 
- ✅ Environment variables
- ✅ Bot token validity

## 🤖 Bot Commands

Once running, your bot supports:

- `/start` - Welcome message
- `/help` - Show all commands
- `/examples` - Chemical compound examples
- `/table` - Periodic table display
- `/ask [element]` - Element information
- `/reactions` - Available animations
- `/animate [reaction]` - Create animation
- `/compounds [search]` - Compound database
- Send any chemical name for conversion

## 🔧 Troubleshooting

### Dependency Issues
```bash
pip install -r requirements_deploy.txt
```

### Token Issues
- Get token from @BotFather
- Set environment variable
- Check for typos

### Import Errors
- Ensure all files are in same directory
- Check Python version (3.8+)
- Verify dependencies installed

## 📦 Deployment Ready

Your bot is ready for deployment on:
- Replit (current environment)
- Heroku
- DigitalOcean
- AWS
- Any Python 3.8+ environment

## 🌟 Features

- **119+ Chemical Compounds** with enhanced formulas
- **37 Reaction Animations** with universal compatibility
- **4K Periodic Table** with detailed element data
- **IUPAC Name Conversion** with 3D visualization
- **Robust Error Handling** with multiple fallbacks

Your Chemistry Telegram Bot is now ready to serve users 24/7!