# Render.com Deployment Guide

## 🚀 Quick Deploy to Render.com

### 1. Repository Setup
1. Upload all files to your GitHub repository
2. Include these Render-specific files:
   - `render.yaml` - Render configuration
   - `runtime.txt` - Python version
   - `Procfile` - Process definition
   - `requirements_deploy.txt` - Dependencies

### 2. Render.com Configuration

1. **Create New Web Service**
   - Connect your GitHub repository
   - Select branch: `main` or `master`
   - Runtime: `Python 3`

2. **Build Settings**
   - Build Command: `pip install -r requirements_deploy.txt`
   - Start Command: `python3 start_bot.py`

3. **Environment Variables**
   - Add `TELEGRAM_BOT_TOKEN` = `your_bot_token_here`
   - Add `LOG_LEVEL` = `INFO`
   - Add `MAX_ANIMATION_SIZE` = `5000000`
   - Add `ANIMATION_TIMEOUT` = `30`

4. **Advanced Settings**
   - Plan: `Starter` (free tier)
   - Region: Choose closest to your users
   - Auto-Deploy: `Yes`

### 3. Deployment Process

1. **Initial Deploy**
   - Click "Create Web Service"
   - Wait for build to complete (5-10 minutes)
   - Check deployment logs

2. **Health Check**
   - Bot will automatically start
   - Check logs for "Bot is ready to receive messages!"
   - Test bot with `/start` command

### 4. Render.yaml Configuration

```yaml
services:
  - type: web
    name: chemistry-telegram-bot
    env: python
    plan: starter
    buildCommand: pip install -r requirements_deploy.txt
    startCommand: python3 start_bot.py
    envVars:
      - key: TELEGRAM_BOT_TOKEN
        sync: false
```

### 5. Dependencies (requirements_deploy.txt)

Optimized for Render.com Python 3.11:
- `requests==2.32.3`
- `python-telegram-bot==21.8`
- `rdkit==2024.3.6`
- `matplotlib==3.8.4`
- `numpy==1.26.4`
- `Pillow==10.4.0`
- `pandas==2.2.2`

### 6. Troubleshooting

**Build Failures:**
- Check Python version in `runtime.txt`
- Verify dependency versions in `requirements_deploy.txt`
- Review build logs for specific errors

**Bot Not Responding:**
- Verify `TELEGRAM_BOT_TOKEN` is set correctly
- Check application logs for errors
- Ensure bot is running with `/start` test

**Memory Issues:**
- Upgrade to paid plan if needed
- Monitor resource usage
- Optimize image generation settings

### 7. Monitoring

- **Logs**: Check Render dashboard for real-time logs
- **Health**: Bot automatically restarts on failure
- **Performance**: Monitor response times and errors

### 8. Scaling

- **Free Tier**: 512MB RAM, sufficient for basic usage
- **Paid Plans**: More RAM for heavy animation generation
- **Auto-scaling**: Render handles traffic spikes automatically

### 9. Custom Domain (Optional)

1. Add custom domain in Render dashboard
2. Configure DNS settings
3. SSL certificates are automatic

### 10. Backup & Recovery

- Code backed up in GitHub
- Environment variables in Render dashboard
- Database-less architecture simplifies recovery

## 🔧 Files Required for Render.com

✅ `render.yaml` - Service configuration
✅ `runtime.txt` - Python version
✅ `Procfile` - Start command
✅ `requirements_deploy.txt` - Dependencies
✅ `start_bot.py` - Enhanced starter with error handling
✅ All bot Python files
✅ `assets/` folder with images

## 📱 Post-Deployment

1. Test all bot commands
2. Verify animations work
3. Check periodic table display
4. Monitor logs for errors

Your Chemistry Telegram Bot will be live at: `https://your-service-name.onrender.com`

The bot runs 24/7 and handles automatic restarts, scaling, and SSL certificates.