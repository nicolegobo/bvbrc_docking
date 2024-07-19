import base64
import glob
import json
import os
import pandas as pd
import sys


# Function to generate the HTML table rows - center aligns numeric values
def generate_table_html(df):
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

    return headers, rows

def image_to_base64(file_path):
    """Read binary file and return base64-encoded string"""
    with open(file_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string

def parse_ligand_file_to_dict(ligand_filename):
    ligand_dict = {}
    with open(ligand_filename, "r") as file:
        for line in file:
            if line.strip():  # check if line is not empty
                key, value = line.split('\t')  # assuming columns are separated by tab
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
    report_subligand_df = pd.DataFrame()
    # Set up variables from the input details dict
    url_base = input_details_dict["url_base"]
    structure_base = input_details_dict["structure_base"]
    workspace_output_path = input_details_dict["params"]["output_path"]
    workspace_output_name = input_details_dict["params"]["output_file"]
    input_pdb = input_details_dict["params"]["input_pdb"][0]
    #### Loop thorugh subligands in the results csv ###
    for sample in input_details_dict["sample_results"]:
        sml_str = ""
        url_list = []
        subligand_df = pd.read_csv(sample, sep="\t")
        new_cols = ["indent", "Rank", "DiffDock confidence", "lig_sdf", "comb_pdb", "CNN Score", "CNN Affinity", "Vinardo"]
        subligand_df.columns = new_cols
        if subligand_df.shape[0] == 0:
            sys.stderr.write("No output for: {} \n".format(sample))
            continue
        else:
        ### START: Edit the rows for each sub ligand dataframes to edit the protein viewer links ###
            for index, row in subligand_df.iterrows():
                #workspace_output_name = input_details_dict["workspace_output_name"].replace(" ", "%20")
                workspace_output_name.replace(" ", "%20")
                comb_pdb = row["comb_pdb"]
                sample_pdb = row["indent"]
                ### START: Make the ligand subtables HTML ###
                # Add protein viewer links for top 3 ranking
                sml_str = input_ligand_dict[sample_pdb]
                viewer_URL = f"{structure_base}={workspace_output_path}/.{workspace_output_name}/{input_pdb}/{sample_pdb}/{comb_pdb}"
                link_html = '<a href="{}" target="_blank">Structure</a>'.format(viewer_URL)
                url_list.append(link_html)
                link_html = ""
            sample_pdb_text = sample_pdb
        # Now back in the sample level for loop #
        subligand_df["Viewer"] = url_list
        sml_str = input_ligand_dict[sample_pdb_text]
        report_subligand_df = subligand_df[["Rank", "Viewer", "Vinardo", "DiffDock confidence", "CNN Score", "CNN Affinity"]]
        # link to subligan directory on workspace
        ligand_dir = f"{url_base}/workspace{workspace_output_path}/.{workspace_output_name}/{input_pdb}/{sample_pdb}"
        ### Write ligand subtables HTML ###
        ligand_subtable_headers, ligand_subtable_rows = generate_table_html(report_subligand_df)


        ligand_subtable_html = \
            """
            <h3 id="{sample_pdb_text}">
                <a href="{ligand_dir}" target="_blank">{sample_pdb_text}</a>: {sml_str}
            </h3>
                <p>Ranked docking conformations:<p>
                    <table>
                        <thead>
                            <tr>
                                {ligand_subtable_headers}
                            </tr>
                        </thead>
                        <tbody>
                            {ligand_subtable_rows}
                        </tbody>
                    </table>
            """.format(sample_pdb_text=sample_pdb_text, \
                       ligand_dir = ligand_dir, \
                        ligand_subtable_headers=ligand_subtable_headers, \
                        ligand_subtable_rows=ligand_subtable_rows, \
                        sml_str = sml_str)

        # add HTML for each table
        ligand_subtables += ligand_subtable_html
        ### END: Make the ligand subtables HTML ###
        ### START: Make the main table HTML ###
        # save the first row for the top ranked ligands table
        first_row = subligand_df.iloc[:1]
        top_ranked_ligands = pd.concat([top_ranked_ligands, first_row], ignore_index=True)

    top_ranked_ligands["SMILES"] = top_ranked_ligands["indent"].map(input_ligand_dict)
    top_ranked_ligands["Ligand ID"] = top_ranked_ligands["indent"].copy()
    top_ranked_ligands["Ligand ID"] = top_ranked_ligands["Ligand ID"].apply(lambda x: f'<a href="#{x}">{x}</a>')
    top_ranked_ligands = top_ranked_ligands[["Ligand ID", "Viewer", "Vinardo", "DiffDock confidence", "CNN Score", "CNN Affinity", "SMILES"]]
    top_ranked_ligands = top_ranked_ligands.sort_values(by="Vinardo")
    ## START: Make the main table HTML ###
    main_ligand_table_headers, main_ligand_table_rows = generate_table_html(top_ranked_ligands)
    # Make the top ligand main table
    main_table_raw_html = \
                """
                    <h3> Successfully Docked Ligands </h3>
                    <p>Following is the top-ranked conformation for each ligand that docked successfully with \
                    the protein. The ligand ID is linked to the list of all conformations for that docking result. \
                    The "structure" link will display a structure viewer of the ligand docked to the protein. <p>
                        <table>
                            <thead>
                                <tr>
                                    {main_ligand_table_headers}
                                </tr>
                            </thead>
                            <tbody>
                                    {main_ligand_table_rows}
                            </tbody>
                        </table>
    """.format(main_ligand_table_headers = main_ligand_table_headers, main_ligand_table_rows = main_ligand_table_rows)
    ### END: Take top row for main table ###
    return {"ligand_subtable_raw_html" : ligand_subtables, "main_table_html" : main_table_raw_html}


def write_html_report(bvbrc_logo_path, main_table, ligand_subtables, input_details_dict):
    base64_string = image_to_base64(bvbrc_logo_path)
    bvbrc_logo_base64 = f'<div class="image-container"><img src="data:image/png;base64,{base64_string}" alt="Embedded Image"></div>'
    protein_title = input_details_dict["proteins"][0]["title"]
    input_pdb = input_details_dict["params"]["input_pdb"][0]

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
        Explore the protein-ligand interactions through advanced artifical intelligence (AL) and machine learning (ML).  This service utilizes a diffusion model,\
            <a href="https://arxiv.org/abs/2210.01776" target="_blank">DiffDock</a>
        to compute a set of poses for a target protein structure and a set of small-molecule ligands. The aim is to simulate and analyze potential binding scenarios \
        “in silico”.  Offering a crucial advantage by predicting the success of protein ligand combinations ahead of costly and time-consuming in vivo experiments. \
        This report summarizes the predicted and evaluated molecular interactions. Clicking structure links out to our protein structure viewer. Clicking the ligand \
        name will jump down that specific ligand in the per-ligand details section of the report. \
        
        <h3>About the Analysis Workflow</h3>
        This service utilizes DiffDock to purpose a set of docking configurations for a given protein structure and \
        a single small molecule.  The diffusion model optimizes the positions and orientations of the small molecule acround the target protein. \
        The sampeled binding poses are then passed to a confidence model which scores and ranks each pose. To compare across different ligands, \
        the binding psoes are then scored with GNINA that contains both the CNN modules and Vinardo score function. \

        In our setup, each ligand in DiffDock's top three docked poses are scored with the GINA moel and ranked based on the Vinardo score.  \
        For a more detailed description of the models, please visit DiffDock, GNINA and Vinardo. \
        <p>

            {main_table}
 <p>The columns are as follows<p>
                <ul>
                    <li><strong>Ligand ID</strong>: The name of the ligand. If a ligand ID is not provided the program will \
                        assign “ID #” where the number is assigned according to the line number in the input.</li>
                    <br>
                    <li><strong>Viewer</strong>: Click on the hyper link to open a new table with the BVBRC protein structure \
                        viewer showing that specific protein ligand interaction.</li>
                    <br>
                    <li><strong>DiffDocK Confidence</strong>: Confidence score of this result as assigned by DiffDock. A lower confidence score \
                         indicates more confidence in the ligand protein docking successfully. </li>
                    <br>
                    <li><strong>CNN Score</strong>: CNN refers to  a type of artificial intelligence called convolutional neural network \
                        (CNNs).  The CNN Score gives a numerical value that represents how how plausible is the binding pose base on the \
                        CNN model evaluation. A higher score indicates better \
                        docking potential. </li>
                    <br>
                    <li><strong>CNN Affinity</strong>: CNN affinity is a *hypothetical* measurement of the ligands affinity to dock with \
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
        main_table=main_table,
        )   
    return html_template

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
    ligand_file = input_details_dict["params"]["ligand_file"]
    input_ligand_dict = parse_ligand_file_to_dict(ligand_file)
    return input_details_dict, input_ligand_dict


def main(argv):
    #Set up the variables and ligand dictionary
    input_details_dict, input_ligand_dict = report_setup(argv)
    ligand_dict = parse_sample_results(input_details_dict, input_ligand_dict)
    ligand_subtables = ligand_dict["ligand_subtable_raw_html"]
    main_table = ligand_dict["main_table_html"]
    bvbrc_logo_path = input_details_dict["bvbrc_logo"]
    html_template= write_html_report(bvbrc_logo_path, main_table, ligand_subtables, input_details_dict)
    output_dir = input_details_dict["output_dir"]
    html_report_path = os.path.join(output_dir, "small_molecule_docking_report.html")
    with open(html_report_path, 'w') as file:
        file.write(html_template)
    sys.stderr.write("Generated HTML report at {}.".format(html_report_path))

if __name__ == "__main__":
    main(sys.argv)
