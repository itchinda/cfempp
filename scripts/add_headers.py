import os

# Configuration de l'entête
HEADER = """//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************
"""

# Extensions cibles
EXTENSIONS = {'.h', '.tpp', '.cpp', '.cc'}
# Dossiers à ignorer (ex: build, .git)
IGNORE_DIRS = {'build', 'build_release', 'build_debug', '.git', '.github', 'venv', 'docs', 'external'}

def add_header_to_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Vérification si l'entête est déjà présent
    if "// --- CFEM++ Library" in content:
        print(f"[-] Skip: {filepath} (Header already exists)")
        return

    # Ajout de l'entête au début du fichier
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(HEADER + "\n" + content)
        print(f"[+] Added: {filepath}")

def main():
    project_root = "."  # Exécuter à la racine du projet
    for root, dirs, files in os.walk(project_root):
        # Filtrer les dossiers à ignorer
        dirs[:] = [d for d in dirs if d not in IGNORE_DIRS]
        
        for file in files:
            if any(file.endswith(ext) for ext in EXTENSIONS):
                filepath = os.path.join(root, file)
                add_header_to_file(filepath)

if __name__ == "__main__":
    main()
