#!/bin/bash

# Chemistry Telegram Bot Deployment Script
# This script helps deploy the bot to various platforms

set -e

echo "🚀 Chemistry Telegram Bot Deployment Script"
echo "==========================================="

# Check if .env file exists
if [ ! -f ".env" ]; then
    echo "⚠️  .env file not found. Creating from template..."
    cp .env.example .env
    echo "📝 Please edit .env file with your bot token before continuing"
    exit 1
fi

# Check if Python 3.8+ is available
python_version=$(python3 --version | cut -d' ' -f2 | cut -d'.' -f1-2)
required_version="3.8"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then
    echo "❌ Python $required_version or higher is required. Found: $python_version"
    exit 1
fi

echo "✅ Python version check passed: $python_version"

# Create virtual environment
echo "🔨 Creating virtual environment..."
python3 -m venv venv
source venv/bin/activate

# Install dependencies
echo "📦 Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Check if TELEGRAM_BOT_TOKEN is set
if [ -z "$TELEGRAM_BOT_TOKEN" ]; then
    echo "⚠️  TELEGRAM_BOT_TOKEN not found in environment"
    echo "Please set your bot token in .env file"
    exit 1
fi

echo "✅ Environment setup complete!"

# Offer deployment options
echo ""
echo "🚀 Deployment Options:"
echo "1. Run locally (development)"
echo "2. Run with systemd (production)"
echo "3. Docker deployment"
echo "4. Exit"

read -p "Choose an option (1-4): " choice

case $choice in
    1)
        echo "🏃 Running bot locally..."
        python3 run.py
        ;;
    2)
        echo "🔧 Setting up systemd service..."
        # Create systemd service file
        sudo tee /etc/systemd/system/chemistry-bot.service > /dev/null <<EOF
[Unit]
Description=Chemistry Education Telegram Bot
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$(pwd)
Environment=PATH=$(pwd)/venv/bin
ExecStart=$(pwd)/venv/bin/python run.py
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF
        
        sudo systemctl daemon-reload
        sudo systemctl enable chemistry-bot.service
        sudo systemctl start chemistry-bot.service
        
        echo "✅ Systemd service created and started"
        echo "📊 Check status: sudo systemctl status chemistry-bot"
        echo "📋 View logs: sudo journalctl -u chemistry-bot -f"
        ;;
    3)
        echo "🐳 Docker deployment not implemented yet"
        echo "Please use option 1 or 2 for now"
        ;;
    4)
        echo "👋 Goodbye!"
        exit 0
        ;;
    *)
        echo "❌ Invalid option"
        exit 1
        ;;
esac