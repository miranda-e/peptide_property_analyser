import os
import subprocess
import shutil
import pandas as pd

def split_fasta(input_filename, output_folder):
     print("Splitting FASTA files...")
     with open(input_filename, 'r') as input_file:
         data = input_file.read().split('>')[1:]
    
     os.makedirs(output_folder, exist_ok=True)

     for index, protein_data in enumerate(data, start=1):
         protein_lines = protein_data.strip().split('\n')
         protein_name = protein_lines[0].split('|')[1]  # Extract the desired string
         protein_sequence = ''.join(protein_lines[1:])
        
         output_filename = os.path.join(output_folder, f'{protein_name}.fasta')  # Use the extracted string as the filename
         with open(output_filename, 'w') as output_file:
             output_file.write(f'>{protein_name}\n{protein_sequence}\n')
     print("FASTA files split.")

def apply_predict_property(fasta_folder, predict_property_path, output_folder):
     print("Applying Predict_Property...")
     os.makedirs(output_folder, exist_ok=True)
    
     for filename in os.listdir(fasta_folder):
         if filename.endswith('.fasta'):
             fasta_path = os.path.join(fasta_folder, filename)
             output_path = os.path.join(output_folder, filename.replace('.fasta', '_PRED'))
            
             subprocess.run([predict_property_path, '-i', fasta_path, '-o', output_path], check=True)
     print("Predict_Property applied.")

def copy_ss8_files(output_folder, target_directory):
     print("Processing SS8 files...")
     os.makedirs(target_directory, exist_ok=True)
    
     for root, _, files in os.walk(output_folder):
         for file in files:
             if file.endswith('.ss8'):
                 file_path = os.path.join(root, file)
                 shutil.copy(file_path, target_directory)

if __name__ == '__main__':
    # Split FASTA files
    input_filename = 'test_input.fasta'
    output_folder = 'test_split_fastas'
    split_fasta(input_filename, output_folder)
    
    # Apply Predict_Property

    fasta_folder = 'test_split_fastas'
    predict_property_path = os.path.expanduser('~/peptide_property_analyser/Predict_Property/Predict_Property.sh')
    output_folder_pp = 'test2_PP_output'
    apply_predict_property(fasta_folder, predict_property_path, output_folder_pp)
    
    # Copy ss8 files
    output_folder_copy = 'test2_PP_output'
    target_directory = os.path.expanduser('~/peptide_property_analyser/test_ss8_files')
    copy_ss8_files(output_folder_copy, target_directory)
    
    # Process ss8 files and create results
    dataframes_directory = target_directory
    file_list = [file for file in os.listdir(dataframes_directory) if file.endswith(".ss8")]

    substrings_df = pd.read_csv('test_peps.csv', skiprows=1, header=None, names=['substring'])

    dataframes = {}  # Initialize dataframes dictionary

    for file in file_list:
        df_name = os.path.splitext(file)[0]
        df = pd.read_csv(os.path.join(dataframes_directory, file), comment='#', delimiter='\s+', header=None, usecols=[0, 1, 2, 3, 4, 5])
        dataframes[df_name] = df

    results_dict = {}

    for substring in substrings_df['substring']:
        substring_results = {}

        for df_name, df in dataframes.items():
            joined_column = ''.join(df[1])
            substring_index = joined_column.find(substring)
            
            if substring_index != -1:
                info_letters = df.loc[substring_index:substring_index + len(substring) - 1, 2].tolist()
                final_string = ''.join(info_letters)
                substring_results[df_name] = final_string
        
        if not substring_results:  # No matches found
            substring_results['No Match'] = "No Match"
        
        results_dict[substring] = substring_results

    # Create the DataFrame from the results_dict
    df_results = pd.DataFrame(columns=['Substring', 'Results', 'pep_LRP_key'])

    for substring, substring_results in results_dict.items():
        for df_name, result in substring_results.items():
            df_results = df_results.append({'Substring': substring, 'Results': result, 'pep_LRP_key': f"{substring}_{df_name}"}, ignore_index=True)

    # Save the DataFrame to JSON and CSV
    df_results.to_json('test2_results_ss8.json', orient='records', indent=4)
    df_results.to_csv('test2_results_ss8.csv', index=False)

    print("Processed SS8 files and created results CSV.")
import pandas as pd

def process_csv(input_csv, output_csv):
    df = pd.read_csv(input_csv)

    # Add "multiple match" column based on duplicates in "Substring" column
    df['multiple match'] = df.duplicated(subset=['Substring'], keep=False)
    df['multiple match'] = df['multiple match'].map({True: 'yes', False: 'no'})

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)

    # Initialize an empty dictionary to store matching Substrings and their corresponding FinalStrings
    matching_strings = {}

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        if row['multiple match'] == 'yes':
            if row['Substring'] not in matching_strings:
                matching_strings[row['Substring']] = [row['Results']]
            else:
                matching_strings[row['Substring']].append(row['Results'])

    # Initialize empty lists to store values for the new columns
    identical_structures = []
    characters_present = {char: [] for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']}
    percentage_by_character = {char: [] for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']}
    majority_characters = []

    # Initialize dictionaries to store character counts in each substring
    char_counts_by_substring = {char: {} for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']}

    # Iterate through each row again and perform the checks
    for index, row in df.iterrows():
        if row['multiple match'] == 'yes':
            if len(matching_strings[row['Substring']]) > 1:
                # Check if all Results are identical
                if all(final_str == matching_strings[row['Substring']][0] for final_str in matching_strings[row['Substring']]):
                    identical_structures.append('yes')
                else:
                    identical_structures.append('no')
            else:
                identical_structures.append('no')
        else:
            identical_structures.append('')

        # Check for the presence of characters in "Results"
        for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']:
            characters_present[char].append('yes' if char in row['Results'] else 'no')

            # Count instances of characters in "Results"
            char_count = row['Results'].count(char)
            char_counts_by_substring[char][row['Substring']] = char_count

        # Calculate percentages of characters compared to total characters
        total_chars = len(row['Results'])
        for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']:
            percentage_by_character[char].append(char_counts_by_substring[char][row['Substring']] / total_chars * 100)

        # Determine the majority character among characters
        max_count = max(char_counts_by_substring[char][row['Substring']] for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I'])
        majority_chars = [char for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I'] if char_counts_by_substring[char][row['Substring']] == max_count]
        majority_characters.append(','.join(majority_chars))

    # Add the new columns to the DataFrame
    df['structures are identical (if multiple)'] = identical_structures
    for char in ['H', 'E', 'L', 'T', 'S', 'G', 'B', 'I']:
        df[f'{char} present'] = characters_present[char]
        df[f'{char} count'] = [char_counts_by_substring[char][row['Substring']] for _, row in df.iterrows()]
        df[f'percentage {char}'] = percentage_by_character[char]
    df['majority character'] = majority_characters

    # Add the "pep_LRP_key" column
    df['pep_LRP_key'] = df['Substring'] + '_' + df['Results']

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)

if __name__ == '__main__':
    input_csv = 'test2_results_ss8.csv'
    output_csv = 'test2_results_ss8_processed.csv'
    process_csv(input_csv, output_csv)
import os
import shutil
import pandas as pd
import json

def load_dataframes(directory):
    """
    Load DataFrame from .acc files in the given directory.
    """
    dataframes = {}
    for file in os.listdir(directory):
        if file.endswith(".acc"):
            df_name = os.path.splitext(file)[0]
            df = pd.read_csv(os.path.join(directory, file), comment='#', delimiter='\s+', header=None, usecols=[0, 1, 2, 3, 4, 5])
            dataframes[df_name] = df
    return dataframes

def process_acc_files(acc_directory, substrings_csv, output_json):
    """
    Process .acc files and create results in JSON format.
    """
    dataframes = load_dataframes(acc_directory)
    
    substrings_df = pd.read_csv(substrings_csv, skiprows=1, header=None, names=['substring'])
    results_dict = {}

    print("Processing RSA files...")
    for substring in substrings_df['substring']:
        substring_results = {}
        for df_name, df in dataframes.items():
            joined_column = ''.join(df[1])
            substring_index = joined_column.find(substring)
            if substring_index != -1:
                info_letters = df.loc[substring_index:substring_index + len(substring) - 1, 2].tolist()
                final_string = ''.join(info_letters)
                substring_results[df_name] = final_string
        if not substring_results:
            substring_results['No Match'] = "No Match"  # Add "No Match" if no match found
        results_dict[substring] = substring_results

    with open(output_json, 'w') as json_file:
        json.dump(results_dict, json_file, indent=4)

def convert_json_to_csv(json_filename, csv_filename):
    """
    Convert JSON results to CSV.
    """
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)
        rows = []
        for substring, substring_results in data.items():
            for df_name, final_string in substring_results.items():
                pep_LRP_key = f"{substring}_{df_name}"  # Create pep_LRP_key
                rows.append([substring, df_name, final_string, pep_LRP_key])

        df = pd.DataFrame(rows, columns=['Substring', 'DataFrame', 'FinalString', 'pep_LRP_key'])
        df.to_csv(csv_filename, index=False)

if __name__ == '__main__':
    # Set paths and filenames
    acc_directory = 'test2_PP_output'  # Adjust the path accordingly
    substrings_csv = 'test_peps.csv'
    output_json = 'test_results_rsa.json'
    output_csv = 'test_results_rsa.csv'

    # Add a step to copy .acc files to 'all_acc' directory
    all_acc_directory = 'test_acc'
    os.makedirs(all_acc_directory, exist_ok=True)

    for root, dirs, files in os.walk(acc_directory):
        for file in files:
            if file.endswith('.acc'):
                src_path = os.path.join(root, file)
                dest_path = os.path.join(all_acc_directory, file)
                shutil.copy(src_path, dest_path)

    # Process acc files and create results in JSON
    process_acc_files(all_acc_directory, substrings_csv, output_json)

    # Convert JSON results to CSV
    convert_json_to_csv(output_json, output_csv)

import pandas as pd

def process_csv(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    
    # Add "multiple match" column based on duplicates in "Substring" column
    df['multiple match'] = df.duplicated(subset=['Substring'], keep=False)
    df['multiple match'] = df['multiple match'].map({True: 'yes', False: 'no'})
    
    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print("Processed RSA files and created results.")
    
    # Initialize an empty dictionary to store matching Substrings and their corresponding FinalStrings
    matching_strings = {}
    
    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        if row['multiple match'] == 'yes':
            if row['Substring'] not in matching_strings:
                matching_strings[row['Substring']] = [row['FinalString']]
            else:
                matching_strings[row['Substring']].append(row['FinalString'])
    
    # Initialize empty lists to store values for the new columns
    identical_structures = []
    b_present = []
    m_present = []
    e_present = []
    majority_structure = []
    percentage_b = []
    percentage_m = []
    percentage_e = []
    tied_characters = []
    
    # Iterate through each row again and perform the checks
    for index, row in df.iterrows():
        if row['multiple match'] == 'yes':
            if len(matching_strings[row['Substring']]) > 1:
                # Check if all FinalStrings are identical
                if all(final_str == matching_strings[row['Substring']][0] for final_str in matching_strings[row['Substring']]):
                    identical_structures.append('yes')
                else:
                    identical_structures.append('no')
            else:
                identical_structures.append('no')
        else:
            identical_structures.append('')
        
        # Check for the presence of characters in "FinalString"
        b_present.append('yes' if 'B' in row['FinalString'] else 'no')
        m_present.append('yes' if 'M' in row['FinalString'] else 'no')
        e_present.append('yes' if 'E' in row['FinalString'] else 'no')
        
        # Count instances of characters in "FinalString"
        b_count = row['FinalString'].count('B')
        m_count = row['FinalString'].count('M')
        e_count = row['FinalString'].count('E')
        
        # Calculate percentages of "B", "M", and "E" characters compared to total characters
        total_chars = len(row['FinalString'])
        percentage_b.append(b_count / total_chars * 100)
        percentage_m.append(m_count / total_chars * 100)
        percentage_e.append(e_count / total_chars * 100)
        
        # Determine the majority structure among "B", "M", and "E" characters
        counts = {'B': b_count, 'M': m_count, 'E': e_count}
        max_count = max(counts.values())
        tied_chars = [char for char, count in counts.items() if count == max_count]
        majority_structure.append(','.join(tied_chars) if len(tied_chars) > 1 else tied_chars[0])
    
    # Add the new columns to the DataFrame
    df['structures are identical (if multiple)'] = identical_structures
    df['B present'] = b_present
    df['M present'] = m_present
    df['E present'] = e_present
    df['majority structure'] = majority_structure
    df['percentage B'] = percentage_b
    df['percentage M'] = percentage_m
    df['percentage E'] = percentage_e
    
    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)

if __name__ == '__main__':
    input_csv = 'test_results_rsa.csv'
    output_csv = 'test_results_processed_rsa.csv'
    process_csv(input_csv, output_csv)
    print("Peptide Property Analyser is complete")
