from pathlib import Path
import subprocess,re
p=Path('docs/developers.toml');s=p.read_text();a=s.index('path = "tests/packaging/test_drawing_notices.py"');b=s.index('\n[[section.note.scope]]',a);block=s[a:b];digest=subprocess.check_output(['git','hash-object','tests/packaging/test_drawing_notices.py'],text=True).strip();block=re.sub('digest = "[a-f0-9]+"','digest = "'+digest+'"',block);s=s[:a]+block+s[b:];p.write_text(s)
