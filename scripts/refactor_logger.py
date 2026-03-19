import os
import re

# Configurations
TARGET_DIRS = ['./libs'] 
FILE_EXTENSIONS = ('.h', '.tpp', '.cpp', '.hpp', '.cxx')

MACRO_MAP = {
    'CFEM_INFO': 'CFEM_INFO',
    'CFEM_LOG': 'CFEM_INFO',
    'CFEM_WARNING': 'CFEM_WARN',
    'CFEM_ERROR': 'CFEM_ERROR'
}

def refactor_content(content):
    for old_macro, new_macro in MACRO_MAP.items():
        # Updated Regex:
        # (?<!\w)         : Word boundary
        # {old_macro}     : The macro name
        # \s*\(           : Opening paren
        # (.*?)           : Non-greedy capture (includes newlines with DOTALL)
        # \)\s*;          : Closing paren followed by optional whitespace and semicolon
        pattern = rf'(?<!\w){old_macro}\s*\((.*?)\)\s*;'
        
        def replacement_logic(match):
            inner_content = match.group(1)
            
            # 1. Replace '+' with '<<'
            # 2. Handle potential multi-line: strip leading/trailing whitespace 
            #    from each line to keep it clean.
            lines = inner_content.split('\n')
            processed_lines = []
            for line in lines:
                # Replace + with << but be careful not to break logic inside quotes
                # Simple replacement is usually safe for logging macros
                line = line.replace('+', '<<')
                processed_lines.append(line)
            
            new_inner = "\n".join(processed_lines)
            
            # Cleanup: If we have multiple '<<' because of the multiline structure
            # we want to ensure it looks like a clean stream.
            return f"{new_macro} << {new_inner.strip()};"

        # re.DOTALL is the magic key here to allow the '.' to catch newlines
        content = re.sub(pattern, replacement_logic, content, flags=re.DOTALL)
    
    return content

def process_files():
    count = 0
    for target_dir in TARGET_DIRS:
        if not os.path.exists(target_dir):
            continue

        for root, _, files in os.walk(target_dir):
            for file in files:
                if file.endswith(FILE_EXTENSIONS):
                    file_path = os.path.join(root, file)
                    
                    with open(file_path, 'r', encoding='utf-8') as f:
                        original_content = f.read()
                    
                    new_content = refactor_content(original_content)
                    
                    if new_content != original_content:
                        with open(file_path, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"Refactored: {file_path}")
                        count += 1
    print(f"\nRefactoring complete. {count} files updated.")

if __name__ == "__main__":
    process_files()
