#!/usr/bin/env python3
import os

# --- Configuration ---
DIRECTORIES_TO_SCAN = ['src', 'include', 'libs', 'tests'] # Ajoute tes dossiers ici
EXTENSIONS = ('.h', '.hpp', '.cpp', '.cc', '.cxx', '.tpp')

COPYRIGHT_HEADER = [
    "//************************************************************************\n",
    "// --- CFEM++ Library\n",
    "// ---\n",
    "// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.\n",
    "// --- ALL RIGHTS RESERVED\n",
    "// ---\n",
    "// --- This software is protected by international copyright laws.\n",
    "// --- Use, distribution, or modification of this software in any form,\n",
    "// --- source or binary, for personal, academic, or non-commercial \n",
    "// --- purposes is permitted free of charge, provided that this \n",
    "// --- copyright notice and this permission notice appear in all \n",
    "// --- copies and supporting documentation.\n",
    "// ---\n",
    "// --- The authors and contributors provide this code \"as is\" without \n",
    "// --- any express or implied warranty. In no event shall they be \n",
    "// --- held liable for any damages arising from the use of this software.\n",
    "//************************************************************************\n"
]

def fix_copyright(lines):
    """Add copyright notice if not already added."""
    # Si le fichier commence déjà par des étoiles ou 'CFEM++', on assume qu'il a un header
    if lines and len(lines) > 5 and "CFEM++" in "".join(lines[:10]):
        return lines, False
    
    # Sinon, on injecte le header avec une ligne vide
    return COPYRIGHT_HEADER + ["\n"] + lines, True

def fix_includes(lines):
    """Sort and create two groups respectively for external library and local files."""
    start_idx = -1
    end_idx = -1
    includes = []

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#include'):
            if start_idx == -1:
                start_idx = i
            end_idx = i
            includes.append(stripped)
        elif not stripped or stripped.startswith('//') or stripped.startswith('/*') or stripped.startswith('*'):
            continue
        elif stripped.startswith('#pragma') or stripped.startswith('#ifndef') or stripped.startswith('#define'):
            continue
        else:
            if start_idx != -1:
                break

    if not includes:
        return lines, False

    sys_includes = set()
    loc_includes = set()

    for inc in includes:
        if '<' in inc and '>' in inc:
            sys_includes.add(inc)
        else:
            loc_includes.add(inc)

    sys_list = sorted(list(sys_includes), key=lambda x: (x.lower(), x))
    loc_list = sorted(list(loc_includes), key=lambda x: (x.lower(), x))

    ideal_block = []
    if sys_list:
        ideal_block.extend(sys_list)
    if sys_list and loc_list:
        ideal_block.append("") 
    if loc_list:
        ideal_block.extend(loc_list)

    ideal_str = "\n".join(ideal_block) + "\n"

    # Extraction du bloc actuel pour comparaison
    current_clean = []
    for line in lines[start_idx:end_idx+1]:
        s = line.strip()
        if s.startswith('#include'):
            current_clean.append(s)
        elif s == '' and current_clean and current_clean[-1] != '':
            current_clean.append('')
            
    while current_clean and current_clean[-1] == '':
        current_clean.pop()
        
    current_str = "\n".join(current_clean) + "\n"

    if ideal_str != current_str:
        # Remplacement
        new_lines = lines[:start_idx] + [line + "\n" for line in ideal_block] + lines[end_idx+1:]
        return new_lines, True

    return lines, False

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    original_lines = lines.copy()
    
    lines, copyright_fixed = fix_copyright(lines)
    lines, includes_fixed = fix_includes(lines)

    if copyright_fixed or includes_fixed:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(lines)
        
        actions = []
        if copyright_fixed: actions.append("Copyright")
        if includes_fixed: actions.append("Includes")
        print(f" Update ({' & '.join(actions)}) : {filepath}")

def main():
    print("Starting the CFEM++ code standardization...\n")
    for root_dir in DIRECTORIES_TO_SCAN:
        if not os.path.exists(root_dir):
            continue
        for root, _, files in os.walk(root_dir):
            for file in files:
                if file.endswith(EXTENSIONS):
                    filepath = os.path.join(root, file)
                    process_file(filepath)
    print("\n Standardization completed ! Track the changes using 'git diff'.")

if __name__ == "__main__":
    main()
