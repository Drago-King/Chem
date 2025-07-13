# Deployment Guide for Chemistry Telegram Bot

## Quick Deployment Steps

### 1. Get Your Bot Token
1. Message @BotFather on Telegram
2. Create a new bot with `/newbot`
3. Copy the bot token

### 2. Deploy to GitHub
1. Create a new GitHub repository
2. Upload all files from this package
3. Set repository secrets:
   - `TELEGRAM_BOT_TOKEN`: Your bot token

### 3. Environment Setup
```bash
# Clone your repository
git clone https://github.com/yourusername/chemistry-telegram-bot.git
cd chemistry-telegram-bot

# Create environment file
cp .env.example .env
# Edit .env with your bot token

# Install dependencies
pip install -r requirements.txt

# Run the bot
python3 run.py
```

## Platform-Specific Deployment

### Replit
1. Create a new Python repl
2. Upload all files
3. Add Secret: `TELEGRAM_BOT_TOKEN`
4. Run: `python3 run.py`

### Heroku
1. Create new app
2. Connect to GitHub repository
3. Set Config Vars: `TELEGRAM_BOT_TOKEN`
4. Deploy from GitHub

### DigitalOcean
1. Create a droplet (Ubuntu 20.04)
2. Install Python 3.8+
3. Clone repository
4. Set environment variables
5. Use systemd service for auto-restart

### AWS EC2
1. Launch EC2 instance
2. Install dependencies
3. Configure security groups
4. Set up systemd service
5. Use CloudWatch for monitoring

## Production Checklist

### Security
- [ ] Bot token stored securely
- [ ] .env file not committed to Git
- [ ] Regular token rotation
- [ ] Firewall configured

### Monitoring
- [ ] Logging configured
- [ ] Health checks enabled
- [ ] Error alerting set up
- [ ] Performance monitoring

### Backup
- [ ] Code backed up to GitHub
- [ ] Configuration documented
- [ ] Deployment process documented
- [ ] Recovery procedures tested

## Environment Variables

```bash
# Required
TELEGRAM_BOT_TOKEN=your_bot_token_here

# Optional
LOG_LEVEL=INFO
MAX_ANIMATION_SIZE=5000000
ANIMATION_TIMEOUT=30
BOT_NAME=@your_bot_name
```

## Troubleshooting

### Bot Not Responding
1. Check bot token validity
2. Verify internet connection
3. Check Telegram API status
4. Review bot logs

### Memory Issues
1. Monitor RAM usage
2. Adjust animation settings
3. Use swap if needed
4. Consider upgrading server

### Performance Issues
1. Enable caching
2. Optimize database queries
3. Use CDN for assets
4. Monitor response times

## Maintenance

### Regular Tasks
- Update dependencies monthly
- Monitor bot usage
- Review error logs
- Test all features

### Updates
- Pull latest changes from Git
- Restart bot service
- Test functionality
- Monitor for issues

## Support

For deployment issues:
- Check GitHub issues
- Review documentation
- Contact support team
- Join community Discord