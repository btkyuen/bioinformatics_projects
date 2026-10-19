Use to generate a containerized environment to run Claude Code.

To enable a secure git access, copy .gitignore_global to your git environment outside of the container, then run the commands:

`git config --global core.excludesfile ~/.gitignore_global`

This should prevent Claude Code from committing keys and other information that should not be exposed.

Within the container, add the following lines to ~/.claude/settings.json:

`{
  "permissions": {
    "deny": [
      "Bash(git push --force*)",
      "Bash(git reset --hard*)",
      "Bash(rm -rf *)"
    ]
  }
}`

This should prevent Claude Code from running git commands that could have negative repercussions.

To limit Claude Code's access to a single repository, use the following steps:
1. Inside the container, generate a dedicated SSH key: `ssh-keygen -t ed25519 -C "claude-container-key" -f ~/.ssh/claude_deploy_key -N ""`
2. Go to your **GitHub repository on the web** -> **Settings** -> **Deploy keys** -> **Add deploy key**.
3. Paste the contents of `~/.ssh/claude_deploy_key.pub` and **check the box for "Allow write access"**.
4. Configure git inside your container workspace to use this specific key by adding this to your project's `.git/config` (or a local SSH config): 
`[core]
    sshCommand = "ssh -i /home/claude-dev/.ssh/claude_deploy_key -o IdentitiesOnly=yes"
`
