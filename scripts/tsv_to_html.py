import pandas as pd
import json
import sys

def tsv_to_html(tsv_file_path, html_output_path, tsv_output_path):
    """
    Convert a TSV file to an interactive HTML table.

    Args:
        tsv_file_path (str): Path to the input TSV file.
        html_output_path (str): Path to save the generated HTML file.
    """
    try:
        # Load the TSV file into a pandas DataFrame
        raw_data = pd.read_csv(tsv_file_path, sep="\t")
        raw_data.rename(
                columns={
                    "structure_link_html": "Viewer", 
                    "score": "DiffDock confidence", 
                    "CNNscore": "CNN Score", 
                    "CNNaffinity": "CNN Affinity", 
                    "smile_string": "SMILES",
                    "drugbank_database_link_html": "Drugbank Generic Name",

                },
                inplace=True
                )
        if "Drugbank Generic Name" in raw_data.columns:
                raw_data = raw_data[[
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
                raw_data = raw_data[[
                    "Ligand ID", 
                    "Viewer", 
                    "Vinardo", 
                    "DiffDock confidence", 
                    "CNN Score", 
                    "CNN Affinity", 
                    "SMILES",
                ]]
        out_tsv = raw_data[[
                    "Ligand ID",  
                    "Vinardo", 
                    "DiffDock confidence", 
                    "CNN Score", 
                    "CNN Affinity", 
                    "SMILES",
        ]]
        out_tsv.to_csv(tsv_output_path, sep="\t", index=False)
        # Convert DataFrame to JSON in 'records' orientation so each entry is treated as a list
        json_data = raw_data.to_json(orient='records', indent=4)

        # HTML template with embedded JSON
        html_template = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Explore All Docking Results</title>
            <!-- DataTables CSS -->
            <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 20px;
                }}
                table {{
                    width: 100%;
                }}
                th input {{
                    width: 100%;
                    box-sizing: border-box;
                }}
            </style>
        </head>
        <body>
            <h1>Interactive Table Viewer</h1>
            <table id="dataTable" class="display" style="width:100%">
                <thead id="tableHead"></thead>
                <tbody id="tableBody"></tbody>
            </table>

            <!-- Embedded JSON Data -->
            <script>
                const tableData = {json_data};
            </script>

            <!-- jQuery -->
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <!-- DataTables -->
            <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
            <script>
                function populateTable(data) {{
                    const tableHead = document.getElementById('tableHead');
                    const tableBody = document.getElementById('tableBody');

                    // Clear previous content
                    tableHead.innerHTML = '';
                    tableBody.innerHTML = '';

                    // Extract headers from the first row
                    const headers = Object.keys(data[0]);

                    // Create header row and filter row
                    const headRow = document.createElement('tr');
                    const filterRow = document.createElement('tr');

                    headers.forEach((header, index) => {{
                        // Header row
                        const th = document.createElement('th');
                        th.textContent = header;
                        headRow.appendChild(th);

                        // Filter row
                        const filterCell = document.createElement('th');
                        const filterInput = document.createElement('input');
                        filterInput.type = 'text'; // Accept both numeric and string filters
                        filterInput.placeholder = `Filter ${{header}}`;
                        filterInput.dataset.column = index;
                        filterInput.addEventListener('keyup', function () {{
                            const columnIndex = parseInt(this.dataset.column);
                            const filterValue = this.value.trim();

                            // Determine the data type of the column (numeric or string)
                            const isNumeric = data.every(row => !isNaN(parseFloat(row[header])));

                            if (filterValue) {{
                                if (isNumeric) {{
                                    // Numeric threshold filtering
                                    const parsedValue = parseFloat(filterValue);
                                    if (!isNaN(parsedValue)) {{
                                        $.fn.dataTable.ext.search = [];
                                        $.fn.dataTable.ext.search.push((settings, row) => {{
                                            const cellValue = parseFloat(row[columnIndex]);
                                            return cellValue >= parsedValue; // Filter values >= threshold
                                        }});
                                    }}
                                }} else {{
                                    // String filtering (case-insensitive substring match)
                                    $.fn.dataTable.ext.search = [];
                                    $.fn.dataTable.ext.search.push((settings, row) => {{
                                        const cellValue = row[columnIndex].toLowerCase();
                                        return cellValue.includes(filterValue.toLowerCase());
                                    }});
                                }}
                            }} else {{
                                // Clear filter if the input is empty
                                $.fn.dataTable.ext.search = [];
                            }}

                            // Redraw the table to apply the new filter
                            $('#dataTable').DataTable().draw();
                        }});
                        filterCell.appendChild(filterInput);
                        filterRow.appendChild(filterCell);
                    }});

                    tableHead.appendChild(headRow);
                    tableHead.appendChild(filterRow);

                    // Create data rows
                    data.forEach(row => {{
                        const tr = document.createElement('tr');
                        headers.forEach(header => {{
                            const td = document.createElement('td');
                            td.innerHTML = row[header]; // Render HTML links directly
                            tr.appendChild(td);
                        }});
                        tableBody.appendChild(tr);
                    }});

                    // Initialize DataTables
                    $('#dataTable').DataTable({{
                        pageLength: 25,  // Rows per page
                        lengthMenu: [10, 25, 50, 100],  // Pagination options
                        orderCellsTop: true,  // Keep filters at the top
                        initComplete: function () {{
                            // Automatically focus on the filter input of the first column for usability
                            $('#dataTable thead input').first().focus();
                        }}
                    }});
                }}

                // Populate the table on page load
                document.addEventListener('DOMContentLoaded', function () {{
                    populateTable(tableData);
                }});
            </script>
        </body>
        </html>
        """.format(json_data=json_data)

        # Save the HTML file
        with open(html_output_path, mode='w') as html_file:
            html_file.write(html_template)

        print("HTML file successfully created at {html_output_path}".format(html_output_path=html_output_path))

    except Exception as e:
        print("Error: {e}".format(e=e))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 tsv_to_html.py <input_tsv_file> <output_html_file> <output_tsv_file>")
    else:
        tsv_file_path = sys.argv[1]
        tsv_df =pd.read_csv(tsv_file_path, sep="\t")
        html_output_path = sys.argv[2]
        tsv_output_path = sys.argv[3]
        tsv_to_html(tsv_file_path, html_output_path, tsv_output_path)

