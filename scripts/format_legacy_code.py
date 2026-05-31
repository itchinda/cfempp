#!/usr/bin/env python3

import os
import sys

DIRECTORIES_TO_SCAN = ['apps', 'libs', 'tests']
EXTENSIONS = ('.h', '.hpp', '.cpp', '.cc', '.cxx', '.tpp')
TEMPLATE_PATH = ".github/COPYRIGHT_TEMPLATE.txt"

def load_copyright_template():
    """Reads the copyright template from the external file."""
    if not os.path.exists(TEMPLATE_PATH):
        print(f"ERROR: Copyright template not found at {TEMPLATE_PATH}")
        sys.exit(1)
    with open(TEMPLATE_PATH, 'r', encoding='utf-8') as f:
        # We read lines and keep the trailing newlines
        return f.readlines()

# Load the template globally
COPYRIGHT_HEADER = load_copyright_template()

def fix_copyright(lines):
    """Checks if the copyright header is present and identical to the template."""
    header_len = len(COPYRIGHT_HEADER)
    
    # If the file is shorter than the header, inject the header
    if len(lines) < header_len:
        return COPYRIGHT_HEADER + ["\n"] + lines, True
    
    # Compare line by line, ignoring trailing whitespaces to be robust
    header_is_correct = True
    for i in range(header_len):
        if lines[i].rstrip() != COPYRIGHT_HEADER[i].rstrip():
            header_is_correct = False
            break
            
    if header_is_correct:
        return lines, False
    
    # If incorrect, clean old header/comments at the top and insert new one
    start_point = 0
    while start_point < len(lines) and (lines[start_point].startswith("//") or lines[start_point].strip() == ""):
        start_point += 1
        
    return COPYRIGHT_HEADER + ["\n"] + lines[start_point:], True

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
