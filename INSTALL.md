# Installation Guide

## Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/chemistry-telegram-bot.git
   cd chemistry-telegram-bot
   ```

2. **Set up environment**
   ```bash
   chmod +x deploy.sh
   ./deploy.sh
   ```

3. **Configure bot token**
   - Copy `.env.example` to `.env`
   - Add your Telegram bot token
   - Run the bot

## Manual Installation

### Prerequisites
- Python 3.8 or higher
- Telegram Bot Token (get from @BotFather)
- 2GB+ RAM (for molecular visualization)
- Internet connection (for chemical data APIs)

### Step-by-Step Setup

1. **Install Python dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Set environment variables**
   ```bash
   cp .env.example .env
   # Edit .env with your bot token
   ```

3. **Test the installation**
   ```bash
   python3 -c "from telegram_bot_simple import SimpleTelegramBot; print('✅ Installation successful')"
   ```

4. **Run the bot**
   ```bash
   python3 telegram_bot_simple.py
   ```

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `TELEGRAM_BOT_TOKEN` | Yes | Your Telegram bot token |
| `LOG_LEVEL` | No | Logging level (INFO, DEBUG, ERROR) |
| `MAX_ANIMATION_SIZE` | No | Maximum animation file size (bytes) |
| `ANIMATION_TIMEOUT` | No | Animation generation timeout (seconds) |

## Deployment Options

### Local Development
```bash
python3 run.py
```

### Production (systemd)
```bash
sudo ./deploy.sh
# Choose option 2
```

### Docker (Coming Soon)
```bash
docker build -t chemistry-bot .
docker run -d --env-file .env chemistry-bot
```

## Troubleshooting

### Common Issues

1. **Import Error: RDKit not found**
   ```bash
   pip install rdkit-pypi
   ```

2. **Animation generation fails**
   - Increase memory allocation
   - Check matplotlib installation
   - Verify PIL/Pillow version

3. **Bot token invalid**
   - Verify token in .env file
   - Check @BotFather for token status
   - Ensure no extra spaces in token

4. **Chemical data not loading**
   - Check internet connection
   - Verify API endpoints are accessible
   - Review PubChem/ChEBI availability

### Performance Optimization

1. **Memory Usage**
   - Monitor RAM usage during animations
   - Adjust animation frame count if needed
   - Use swap space if memory is limited

2. **Response Time**
   - Enable local compound database
   - Cache frequently requested structures
   - Use CDN for static assets

## System Requirements

### Minimum
- Python 3.8+
- 1GB RAM
- 500MB storage
- Internet connection

### Recommended
- Python 3.10+
- 2GB+ RAM
- 1GB storage
- Stable internet connection

## Platform Support

### Tested Platforms
- ✅ Linux (Ubuntu 20.04+)
- ✅ macOS (10.15+)
- ✅ Windows 10+
- ✅ Replit
- ✅ Heroku
- ✅ DigitalOcean

### Cloud Deployment
- AWS EC2
- Google Cloud Platform
- Microsoft Azure
- DigitalOcean Droplets
- Heroku (with buildpacks)

## Security Considerations

1. **Bot Token Security**
   - Never commit .env files
   - Use environment variables in production
   - Rotate tokens regularly

2. **API Rate Limits**
   - Respect PubChem API limits
   - Implement request queuing
   - Use caching for frequent requests

3. **User Data**
   - No user data is stored
   - All interactions are stateless
   - No personal information collected

## Support

For installation issues:
1. Check the troubleshooting section
2. Search existing GitHub issues
3. Create a new issue with details
4. Join our Discord server (link in README)

## License

MIT License - see LICENSE file for details