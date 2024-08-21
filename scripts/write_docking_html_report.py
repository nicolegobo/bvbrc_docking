import base64
import glob
import json
import os
import pandas as pd
import sys

from pathlib import Path


def check_RDKit_invalid_ligands(input_details_dict):
    lig_failed_rdkit = os.path.join(input_details_dict["staging_dir"], "invalid_smile_strings.txt")
    # assuming check_input_smile_strings writes out an empty file
    path = Path(lig_failed_rdkit)
    if path.is_file() and path.stat().st_size > 0:
        invalid_ligs_df = pd.read_csv(lig_failed_rdkit, sep='\t', header=None, names=['Ligand ID', 'SMILE'])
        lig_failed_rdkit_table_html = generate_table_html_2(invalid_ligs_df, table_width='95%')
        rdkit_failed_ligands_html = """
        <h3> RDKit Non-compliant Ligands </h3>
        <p>The following ligands did not pass SMILE string validation preformed by <a href="https://rdkit.org/">RDKit</a>. \
            This is a collection of cheminformatics and machine-learning software. RDKit validates SMILE strings according to \
            parsing and sanitization.  For parsing, RDKit uses grammar defined in the <a href="https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/SmilesParse">Smile Parse</a> \
            module. Which closely follows guidelines established in <a href="http://opensmiles.org/opensmiles.html">OpenSmiles</a>. \
            A detailed explaination of sanitization is available <a href="https://www.rdkit.org/docs/RDKit_Book.html">RDKit Book</a> \
            under the Molecular Sanitization header. <p>
            {lig_failed_rdkit_table_html}
        """.format(
        lig_failed_rdkit_table_html=lig_failed_rdkit_table_html,
        )
    else:
        rdkit_failed_ligands_html = ""
    return rdkit_failed_ligands_html

def check_dd_invalid_ligands(input_details_dict, input_ligand_dict):
    # check if the bad ligands file exists
    lig_failed_diffdock = os.path.join(input_details_dict["work_dir"], input_details_dict["params"]["input_pdb"][0], "out", "bad-ligands.txt")
    # add an if file exists because this function may be used if all ligands are invalid and diff dock never runs
    path = Path(lig_failed_diffdock)
    print(path)
    # Check if the file exists and is not empty
    if path.is_file() and path.stat().st_size > 0:
        print("DD failure")
        dd_failed_ligands = []
        with open(lig_failed_diffdock) as file:
            for line in file:
                if line.startswith('@BodyException'):
                    ligand_id = line.split()[1]
                    smile = input_ligand_dict[ligand_id]
                    # Append the key-value pair as a dictionary to the list
                    dd_failed_ligands.append({'Ligand ID': ligand_id, 'Smile': smile})
                    # match with the smile string to match the other tables
        dd_failed_ligands_df = pd.DataFrame(dd_failed_ligands)
        lig_failed_DD_table_html = generate_table_html_2(dd_failed_ligands_df, table_width='95%')
        dd_failed_ligands_html = """
        <h3> DiffDock Undocked Ligands </h3>
        <p> Ligands in this table did not dock to the protein. This coud be because they are incompatible with the protein or the current version of DiffDock. \
        Another reason could be the available memory during your job. To test this, please submit a new job with each ligand invidiaully or in smaller groups. \
        These ligands are also described in the file "bad-ligands.txt" in the protein subdirectory for your job.
        <br> <br>
        If you have questions about failed ligands, we encourage you to reach out to a team member by either reporting the job or contacting us by clicking "About" \
        in our header then the dropdown option "Contact Us".<p>
        {lig_failed_DD_table_html}
        """.format(
        lig_failed_DD_table_html=lig_failed_DD_table_html,
        )
    else:
        dd_failed_ligands_html = ""
    return dd_failed_ligands_html

# Function to generate the HTML table with a standard width
def generate_table_html_2(df, table_width='95%'):
    # Generate table headers
    headers = ''.join(f'<th>{header}</th>' for header in df.columns)

    # Function to determine if a value is numeric
    def is_numeric(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    # Initialize rows string
    rows = ''
    
    # Generate table rows
    for _, row in df.iterrows():
        row_html = ''
        for column in df.columns:
            cell_value = row[column]
            if is_numeric(cell_value):
                row_html += f'<td style="text-align: center;">{cell_value}</td>'
            else:
                row_html += f'<td>{cell_value}</td>'
        rows += f'<tr>{row_html}</tr>'

    # Construct the complete HTML table with specified width
    table_html = f'''
    <table style="width: {table_width}; border-collapse: collapse; border: 1px solid black;">
        <thead>
            <tr>{headers}</tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>
    '''
    return table_html

def get_name_by_id(data_dict, id):
    return data_dict.get(id, "ID not found")

def image_to_base64(file_path):
    """Read binary file and return base64-encoded string"""
    with open(file_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string

def load_data_to_dict(file_path):
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) == 3:
                ligand_id, name, smile_string = columns
                data_dict[ligand_id] = name
    return data_dict

def parse_ligand_file_to_dict(ligand_filename):
    ligand_dict = {}
    with open(ligand_filename, "r") as file:
        for line in file:
            if line.strip():  # check if line is not empty
                key, value = line.split('\t')  # Smile ligand file columns are separated by tab
                ligand_dict[key.strip()] = value.strip()
    return ligand_dict

def parse_ligand_list_to_dict(ligand_smiles_list):
    ligand_dict = {}
    for item in ligand_smiles_list:
        if len(item) == 2:
            key = item[0]
            value = item[1]
            ligand_dict[key] = value
    return ligand_dict

# Each sub directory has a sample level 'results csv' that contains information for each ligand
def parse_sample_results(input_details_dict, input_ligand_dict):
    # All the subtable html will be apended to this string 
    ligand_subtables = ""
    # Set up dataframes to be written to HTML table
    top_ranked_ligands = pd.DataFrame()
    # report_subligand_df = pd.DataFrame()
    # Set up variables from the input details dict
    url_base = input_details_dict["url_base"]
    structure_base = input_details_dict["structure_base"]
    workspace_output_path = input_details_dict["params"]["output_path"]
    workspace_output_name = input_details_dict["params"]["output_file"]
    workspace_output_name_url = workspace_output_name.replace(" ", "%20")
    data = input_details_dict["results"]
    # Transform data into a dataframe
    rows = []
    for protein_id, level2_data in data.items():
        for ligand, records in level2_data.items():
            for record in records:
                row = {'protein': protein_id, "Ligand ID": ligand}
                row.update(record)
                rows.append(row)

    if len(rows) == 0:
        write_html_report_all_ligands_invalid(input_details_dict, input_ligand_dict)
        sys.exit(0)
    dff = pd.DataFrame(rows)
    dff = dff.astype({
                    "Vinardo": float, 
                     "score": float, 
                     "CNNscore": float, 
                     "CNNaffinity": float, 
                    }
                     )
    dff = dff.round(3)
    if input_details_dict["params"]["ligand_library_type"] == "named_library":
        ligand_lib_info = os.path.join(input_details_dict["staging_dir"], "info.txt")
        ligand_lib_dict = load_data_to_dict(ligand_lib_info)
        # Adding the 'names' column based on the ligand id
        dff['Names'] = dff["Ligand ID"].map(ligand_lib_dict)
        dff["drugbank_database_URL"] = "https://go.drugbank.com/drugs/" + dff["Ligand ID"]
        dff["drugbank_database_link_html"] ='<a href="' + dff["drugbank_database_URL"] + '" target="_blank">' + dff['Names'] + '</a>'
    dff["smile_string"] = dff["Ligand ID"].map(input_ligand_dict)
    # path to ligand directory on workspace
    dff["ligand_dir_path"] = f"{url_base}/workspace{workspace_output_path}/.{workspace_output_name_url}/" + dff["protein"] + "/" + dff["Ligand ID"]
    dff["viewer_URL"] = "{}={}/.{}/".format(structure_base, workspace_output_path,workspace_output_name_url) + dff["protein"] + "/" + dff["Ligand ID"] + "/" + dff["comb_pdb"]
    dff["structure_link_html"] = '<a href="' + dff["viewer_URL"] + '" target="_blank">Structure</a>'
    ligand_subtables = ""
    for ligand in pd.unique(dff["Ligand ID"]):
        report_tmp = dff.loc[dff["Ligand ID"] == ligand].copy()
        ligand_id_text = report_tmp["Ligand ID"].iloc[0]
        ligand_dir = report_tmp["ligand_dir_path"].iloc[0]
        sml_str = report_tmp["smile_string"].iloc[0]
        # rename the columns
        report_tmp.rename(
            columns={
                "structure_link_html": "Viewer", 
                "score": "DiffDock confidence", 
                "CNNscore": "CNN Score", 
                "CNNaffinity": "CNN Affinity", 
                "smile_string": "SMILES",
                "drugbank_database_link_html": "Ligand Name",
            },
            inplace=True
        )
        report_tmp = report_tmp.sort_values(by="Vinardo", ascending=True)
        # if "Ligand Name" in report_tmp.columns:
        if "Drugbank Generic Name" in report_tmp.columns:
            report_tmp = report_tmp[[
                "Ligand ID",  
                "Drugbank Generic Name",
                "Viewer", 
                "Vinardo", 
                "DiffDock confidence", 
                "CNN Score", 
                "CNN Affinity", 
                "SMILES",
            ]]
        else:
            report_tmp = report_tmp[[
                "Ligand ID", 
                "Viewer", 
                "Vinardo", 
                "DiffDock confidence", 
                "CNN Score", 
                "CNN Affinity", 
                "SMILES",
            ]]
        one_ligand_subtable_html = generate_table_html_2(report_tmp, table_width='95%')
        ligand_subtable_html = \
            """
            <h3 id="{ligand_id_text}">
                <a href="{ligand_dir}" target="_blank">{ligand_id_text}</a>: {sml_str}
            </h3>
                <p>Ranked docking conformations:<p>
                    {one_ligand_subtable_html}
            """.format(ligand_id_text=ligand_id_text, \
                    ligand_dir = ligand_dir, \
                    one_ligand_subtable_html=one_ligand_subtable_html, \
                        # ligand_subtable_headers=ligand_subtable_headers, \
                        # ligand_subtable_rows=ligand_subtable_rows, \
                        sml_str = sml_str)
        ligand_subtables += ligand_subtable_html
        ### END: Make the ligand subtables HTML ###
        #### START: Make the main table HTML ###
        # save the first row for the top ranked ligands table
        top_ranked_ligands = pd.concat([top_ranked_ligands, report_tmp.iloc[:1]], ignore_index=True)
    report_top_ranked_ligands = top_ranked_ligands.copy()
    report_top_ranked_ligands["Ligand ID"] = '<a href="#' + top_ranked_ligands["Ligand ID"] + '">' + top_ranked_ligands["Ligand ID"]+ '</a>'
    report_top_ranked_ligands = report_top_ranked_ligands.sort_values(by="Vinardo", ascending=True)
    ## START: Make the main table HTML ###
    main_ligand_table_html = generate_table_html_2(report_top_ranked_ligands, table_width='95%')
    # Make the top ligand main table
    main_table_raw_html = \
                """
                    <h3> Successfully Docked Ligands </h3>
                    <p>Following is the top-ranked conformation for each ligand that docked successfully with \
                    the protein. The ligand ID is linked to the list of all conformations for that docking result. \
                    The "structure" link will display a structure viewer of the ligand docked to the protein. <p>
                    {main_ligand_table_html}
    """.format(main_ligand_table_html=main_ligand_table_html)
    ### END: Take top row for main table ###
    return {"ligand_subtable_raw_html" : ligand_subtables, "main_table_html" : main_table_raw_html}

def write_html_report(bvbrc_logo_path, main_table, ligand_subtables, input_details_dict, input_ligand_dict):
    base64_string = image_to_base64(bvbrc_logo_path)
    bvbrc_logo_base64 = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
    protein_title = input_details_dict["proteins"][0]["title"]
    input_pdb = input_details_dict["params"]["input_pdb"][0]
    # check for invalid ligands
    rdkit_failed_ligands_html = check_RDKit_invalid_ligands(input_details_dict)
    dd_failed_ligands_html = check_dd_invalid_ligands(input_details_dict, input_ligand_dict)
    ### write the report HTML ###
    html_template = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <style>
                body {{ font-family: Roboto, sans-serif; color: black; }}
                    header {{
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        padding: 10px 20px;
                    }}
                    header > a img {{
                        max-width: 225px;  /* Maximum width */
                        max-height: 225px;  /* Maximum height */
                        width: auto;
                        height: auto;
                    }}
                    .title {{
                        font-size: 36px;  /* Adjust the size of the title text */
                        font-family: 'Roboto', sans-serif;
                        font-weight: bold;
                        color: black;
                    }}
                    .warning {{ color: black; }}
                    table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
                    th, td {{ padding: 5px; text-align: left; }}
                    img {{ width: 100%; max-width: 600px; height: auto; }}
                    .image-row {{
                        display: flex;
                        flex-wrap: wrap;
                        justify-content: flex-start;
                    }}
                    .image-container {{
                        width: 33%; /* Each image container takes up one-third of the row */
                        padding: 5px; /* Padding around the images */
                        box-sizing: border-box;
                    }}
                    .img {{
                        width: 100%; /* Make the image expand to fill the container */
                        max-width: 600px; /* Maximum width of the image */
            </style>
        </head>
        <body>
            <header>
                <div class="title">Small-molecule Docking Service Report</div>
                <a href="https://www.bv-brc.org/" target="_blank">
                    {bvbrc_logo_base64}
                    </a>
            </header>
        <p> 
        Explore protein-ligand interactions through advanced artificial intelligence (AI) and machine learning (ML).  This service utilizes a diffusion model,\
            <a href="https://arxiv.org/abs/2210.01776" target="_blank">DiffDock</a>
        to compute a set of poses for a target protein structure and a set of small-molecule ligands. The aim is to simulate and analyze potential binding scenarios \
        “in silico”.  Offering a crucial advantage by predicting the success of protein ligand combinations ahead of costly and time-consuming in vivo experiments. \
        This report summarizes the predicted and evaluated molecular interactions. Clicking structure links out to our protein structure viewer. Clicking the ligand \
        name will jump down that specific ligand in the per-ligand details section of the report. \
        
        <h3>About the Analysis Workflow</h3>
        This service utilizes DiffDock to purpose a set of docking configurations for a given protein structure and \
        a single small molecule.  The diffusion model optimizes the positions and orientations of the small molecule around the target protein. \
        The sampled binding poses are then passed to a confidence model which scores and ranks each pose. To compare across different ligands, \
        the binding poses are then scored with GNINA that contains both the CNN modules and Vinardo score function. \

        In our setup, each ligand in DiffDock's top three docked poses are scored with the GINA model and ranked based on the Vinardo score.  \
        For a more detailed description of the models, please visit DiffDock, GNINA and Vinardo. \
        <p>

            {main_table}
 <p>The columns are as follows<p>
                <ul>
                    <li><strong>Ligand ID</strong>: The name of the ligand. If a ligand ID is not provided the program will \
                        assign “ID #” where the number is assigned according to the line number in the input.</li>
                    <br>
                    <li><strong>Viewer</strong>: Click on the hyperlink to open a new table with the BVBRC protein structure \
                        viewer showing that specific protein ligand interaction.</li>
                    <br>
                    <li><strong>DiffDock Confidence</strong>: Confidence score of this result as assigned by DiffDock. A lower confidence score \
                         indicates more confidence in the ligand protein docking successfully. </li>
                    <br>
                    <li><strong>CNN Score</strong>: CNN refers to  a type of artificial intelligence called convolutional neural network \
                        (CNNs).  The CNN Score gives a numerical value that represents how how plausible is the binding pose base on the \
                        CNN model evaluation. A higher score indicates better \
                        docking potential. </li>
                    <br>
                    <li><strong>CNN Affinity</strong>: CNN affinity is a hypothetical measurement of the ligands affinity to dock with \
                        the target protein as calculated by the central neural network described above. A higher affinity value indicates a \
                            a better chance of successful ligand docking.</li>
                    <br>
                    <li><strong>Vinardo</strong>: An empirical score function that evaluate the binding pose with terms from Gaussian steric \
                        attractions, quadratic steric repulsions, Lennard-Jones potentials, electrostatic interactions, hydrophobic interactions, \
                        non-hydrophobic interactions, and non-directional hydrogen bonds.  A lower Vinardo score indicates a better chance of ligand \
                        binding. </li>
                    <br>
                    <li><strong>SMILES</strong>: SMILES is the “Simplified Molecular Input Line Entry System,” which is used to translate a \
                        chemical's three-dimensional structure into a string of symbols that is easily understood by computer software.   For \
                        a user friendly explanation please visit this
                         <a href="https://www.epa.gov/sites/default/files/2015-05/documents/appendf.pdf" target="_blank"> link</a>.
                           </li>
                    <br>
                </ul>
            
            <h1> {input_pdb}: {protein_title} </h1>
            <h3> Per-ligand Details </h3>
            {ligand_subtables}

            {rdkit_failed_ligands_html}
            <br>
            {dd_failed_ligands_html}
            <h3>References</h3>
            <ol>
            <li>Olson RD, Assaf R, Brettin T, Conrad N, Cucinell C, Davis JJ, Dempsey DM, Dickerman A, Dietrich EM, Kenyon RW, Kuscuoglu \
                M, Lefkowitz EJ, Lu J, Machi D, Macken C, Mao C, Niewiadomska A, Nguyen M, Olsen GJ, Overbeek JC, Parrello B, Parrello V, \
                Porter JS, Pusch GD, Shukla M, Singh I, Stewart L, Tan G, Thomas C, VanOeffelen M, Vonstein V, Wallace ZS, Warren AS, \
                Wattam AR, Xia F, Yoo H, Zhang Y, Zmasek CM, Scheuermann RH, Stevens RL. Nucleic Acids Res. 2022 Nov 9:gkac1003. \
                doi: 10.1093/nar/gkac1003. Epub ahead of print. PMID: 36350631</a></li>
        
            <li>Corso, Gabriele, Arthur Deng, Benjamin Fry, Nicholas Polizzi, Regina Barzilay, and Tommi Jaakkola. "Deep Confident Steps to New Pockets: Strategies for Docking Generalization." arXiv preprint arXiv:2402.18396 (2024).</a></li>
        
            <li>McNutt, Andrew T., Paul Francoeur, Rishal Aggarwal, Tomohide Masuda, Rocco Meli, Matthew Ragoza, Jocelyn Sunseri, and David Ryan Koes. "GNINA 1.0: molecular docking with deep learning." Journal of cheminformatics 13, no. 1 (2021): 43.</a></li>
            
            <li>Quiroga R, Villarreal MA (2016) Vinardo: A Scoring Function Based on Autodock Vina Improves Scoring, Docking, and Virtual Screening. PLOS ONE 11(5): e0155183. https://doi.org/10.1371/journal.pone.0155183.</a></li>

            <li>Bento, A.P., Hersey, A., Félix, E. et al. An open source chemical structure curation pipeline using RDKit. J Cheminform 12, 51 (2020). https://doi.org/10.1186/s13321-020-00456-1</a></li>

            <li>David Sehnal, Sebastian Bittrich, Mandar Deshpande, Radka Svobodová, Karel Berka, Václav Bazgier, Sameer Velankar, Stephen K Burley, Jaroslav Koča, Alexander S Rose: Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures, Nucleic Acids Research, 2021; https://doi.org/10.1093/nar/gkab314.</a></li>

            <ol>
        <script>
            document.querySelectorAll('a[target="_blank"]').forEach(link => {{
            link.addEventListener('click', function(event) {{
            event.preventDefault();
            window.open(this.href, '_blank');
            }});
            }});
        </script>
    """.format(
        bvbrc_logo_base64=bvbrc_logo_base64,
        ligand_subtables = ligand_subtables,
        input_pdb = input_pdb,
        protein_title = protein_title,
        main_table = main_table,
        rdkit_failed_ligands_html = rdkit_failed_ligands_html,
        dd_failed_ligands_html = dd_failed_ligands_html,
        )
    output_dir = input_details_dict["output_dir"]
    html_report_path = os.path.join(output_dir, "small_molecule_docking_report.html")
    with open(html_report_path, 'w') as file:
        file.write(html_template)
    sys.stderr.write("Generated HTML report at {}.".format(html_report_path))

def write_html_report_all_ligands_invalid(input_details_dict, input_ligand_dict):
# if no ligands successfully bound, write out to an empty report.
    bvbrc_logo_path = input_details_dict["bvbrc_logo"]
    base64_string = image_to_base64(bvbrc_logo_path)
    bvbrc_logo_base64 = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
    rdkit_failed_ligands_html = check_RDKit_invalid_ligands(input_details_dict)
    dd_failed_ligands_html = check_dd_invalid_ligands(input_details_dict, input_ligand_dict)

    if "proteins" in input_details_dict.keys():
        protein_title = input_details_dict["proteins"][0]["title"]
        input_pdb = input_details_dict["params"]["input_pdb"][0]
        protein_text = """
    <p> Zero ligands successfully bound to the given protein:
        <br>
        <br>
    Input Protein: {protein_title}
        <br>
    Input PBD ID: {input_pdb}
        <br>
        """
    else:
        protein_title = ""
        input_pdb = ""
        protein_text= ""
    ## TO DO: in report refactor make this template less redundant 
    html_template_failed_job = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <style>
            body {{ font-family: Roboto, sans-serif; color: black; }}
                header {{
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    padding: 10px 20px;
                }}
                header > a img {{
                    max-width: 225px;  /* Maximum width */
                    max-height: 225px;  /* Maximum height */
                    width: auto;
                    height: auto;
                }}
                .title {{
                    font-size: 36px;  /* Adjust the size of the title text */
                    font-family: 'Roboto', sans-serif;
                    font-weight: bold;
                    color: black;
                }}
                .warning {{ color: black; }}
                table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
                th, td {{ padding: 5px; text-align: left; }}
                img {{ width: 100%; max-width: 600px; height: auto; }}
                .image-row {{
                    display: flex;
                    flex-wrap: wrap;
                    justify-content: flex-start;
                }}
                .image-container {{
                    width: 33%; /* Each image container takes up one-third of the row */
                    padding: 5px; /* Padding around the images */
                    box-sizing: border-box;
                }}
                .img {{
                    width: 100%; /* Make the image expand to fill the container */
                    max-width: 600px; /* Maximum width of the image */
        </style>
    </head>
    <body>
        <header>
            <div class="title">Small-molecule Docking Service Report</div>
            <a href="https://www.bv-brc.org/" target="_blank">
                {bvbrc_logo_base64}
                </a>
        </header>
        {protein_text}
        {rdkit_failed_ligands_html}
        <br>
        {dd_failed_ligands_html}
        <br>
    """.format(
    bvbrc_logo_base64=bvbrc_logo_base64,
    protein_text=protein_text, 
    protein_title=protein_title,
    input_pdb=input_pdb,
    rdkit_failed_ligands_html=rdkit_failed_ligands_html,
    dd_failed_ligands_html=dd_failed_ligands_html,
    )
    output_dir = input_details_dict["output_dir"]
    html_report_path = os.path.join(output_dir, "small_molecule_docking_report.html")
    with open(html_report_path, 'w') as file:
        file.write(html_template_failed_job)
    sys.stderr.write("Generated HTML report at {}. \n".format(html_report_path))
    sys.stderr.write("Zero ligands bound to the protein. \n")

def report_setup(argv):
    ### Parse input data for the report ###
    # Details from Analysis 
    report_data_file = argv[1]
    input_details_dict = {}
    input_ligand_dict = {}
    # Open and read the JSON file
    with open(report_data_file, 'r') as file:
        file_content = file.read()
    # Parse the JSON content
    input_details_dict = json.loads(file_content)
    work_dir = input_details_dict["work_dir"]
    input_pdb = input_details_dict["params"]["input_pdb"][0]
    input_details_dict["sample_results"]  = glob.glob("{}/{}/out/*/result.csv".format(work_dir, input_pdb))
    if len(input_details_dict["sample_results"]) == 0:
        print('All ligands were marked as invalid by check input smile strings')
        ligand_file = input_details_dict["failed_validation"]
        input_ligand_dict = parse_ligand_file_to_dict(ligand_file)
        return input_details_dict, input_ligand_dict
    else:
        ligand_file = input_details_dict["params"]["ligand_file"]
        input_ligand_dict = parse_ligand_file_to_dict(ligand_file)
        return input_details_dict, input_ligand_dict


def main(argv):
    # Set up the variables and ligand dictionary
    input_details_dict, input_ligand_dict = report_setup(argv)
    if len(input_details_dict["sample_results"]) == 0:
        write_html_report_all_ligands_invalid(input_details_dict, input_ligand_dict)
    else:
        ligand_dict = parse_sample_results(input_details_dict, input_ligand_dict)
        ligand_subtables = ligand_dict["ligand_subtable_raw_html"]
        main_table = ligand_dict["main_table_html"]
        bvbrc_logo_path = input_details_dict["bvbrc_logo"]
        write_html_report(bvbrc_logo_path, main_table, ligand_subtables, input_details_dict, input_ligand_dict)

if __name__ == "__main__":
    main(sys.argv)
