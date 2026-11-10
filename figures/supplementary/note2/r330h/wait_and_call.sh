cd "$(dirname "$0")"
while pgrep -u "$USER" -f "gdc-client download -m r330h/manifest.tsv" >/dev/null; do sleep 60; done
tail -5 logs/dl_login.log
bash call_r330h.sh
