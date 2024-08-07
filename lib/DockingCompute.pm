package DockingCompute;

use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Path qw(rmtree make_path);
use File::Spec;
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON;
use Text::CSV qw(csv);
use Getopt::Long::Descriptive;
use P3DataAPI;
use YAML;
use Template;
use Module::Metadata;
use warnings;
use strict;

use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(api smiles_list debug work_dir staging_dir output_dir app params
			     ligand_name));

sub new
{
    my($class) = @_;
    my $self = {};
    bless $self, $class;

    $self->{api} = P3DataAPI->new();
    
    return $self;
}

sub run
{
    my($self, $app, $app_def, $raw_params, $params) = @_;

    #
    # Set up work and staging directories.
    # 
    #

    $self->app($app);
    $self->params($params);

    my $tmp_dir = $params->{_tmpdir};
    $tmp_dir //= getcwd . "/tmp.$$";
    
    my $work_dir = "$tmp_dir/work";
    my $staging_dir = "$tmp_dir/staging";
    my $output_dir = "$tmp_dir/output";

    $self->{work_dir} = $work_dir;
    $self->{staging_dir} = $staging_dir;
    $self->{output_dir} = $output_dir;

    make_path($work_dir, $staging_dir, $output_dir);

    #
    # Find the lingands
    #
    my $ligand_file;
    if ($params->{ligand_library_type} eq 'named_library')
    {
        if ($params->{ligand_named_library} eq 'approved-drugs')
        {
	    $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/drugbank_approved.txt");
        }
        elsif ($params->{ligand_named_library} eq 'experimental_drugs')
        {
	    $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/drugbank_exp_inv.txt");
        }
        elsif ($params->{ligand_named_library} eq 'test')
        {
        $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/test.txt");
        }
        else
        {
        die "Unknown ligand library type selected $params->{ligand_library_type}";
        }

    }
    elsif ($params->{ligand_library_type} eq 'smiles_list')
    {
	$ligand_file = $self->load_ligand_smiles($params->{ligand_smiles_list});
    }
    elsif ($params->{ligand_library_type} eq 'ws_file')
    {
	$ligand_file = $self->load_ligand_ws_file($params->{ligand_ws_file});
    }

    #
    # Add ligand and parameters to self
    #
    $params->{ligand_file} = $ligand_file; # all ligand inputs are written to file
    $self -> {params} = $params;

    #
    # Stage the PDB data
    #

    my @pdb_list = $self->stage_pdb($params->{input_pdb});

    #
    # And compute
    #

    for my $pdb (@pdb_list)
    {
	my $work = "$work_dir/$pdb->{pdb_id}";
	make_path($work);
	chdir($work);
	$self->compute_pdb($pdb, $ligand_file, $work);
    }

    #
    # Write report
    #

    $self->write_report(\@pdb_list);

    $self->save_output_files($app, $output_dir);
    
    if (!$params->{enable_debug})
    {
	# rmtree($tmp_dir);
    }
}

sub write_report
{
    my($self, $pdbs) = @_;
    
    my $url_base = $ENV{P3_BASE_URL} // "https://www.bv-brc.org";
    my %vars = (
        proteins => $pdbs,
		work_dir => $self->{work_dir},
        staging_dir => $self->{staging_dir},
        output_dir => $self->{output_dir},
		ligands => $self->{ligand_name},
		ligand_info => $self->{ligand_info},
		results => $self->{result_data},
		params => $self->params,
		output_folder => $self->params->{output_path} . "/." . $self->params->{output_file},
		url_base => $url_base,
		feature_base => "$url_base/view/Feature",
		structure_base => "$url_base/view/ProteinStructure#path",
        bvbrc_logo => "/vol/bvbrc/production/application-backend/bvbrc_docking/bv-brc-header-logo-bg.png"
			);

	# Convert the hash to a JSON string 
	my $json_text = to_json(\%vars, { pretty => 1 });
	#Define the path to the report_data.json file
	my $report_data_path = File::Spec->catfile($self->{work_dir}, "report_data.json");

	# Write the JSON string to the file
	open(my $fh, '>', $report_data_path) or die "Could not open file '$report_data_path': $!";
	print $fh $json_text;
	close($fh);
	 print "Analysis data written to $report_data_path\n";

	my @cmd = (
		"write_docking_html_report",
		"$report_data_path"
	);

    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Report command failed $?: @cmd";
    }
}

#
# Compute one PDB.
#
# We write output for this PDB into $output_dir/$pdb_id.
#
# We've already chdir'd to $work_dir.
#

sub compute_pdb
{
    my($self, $pdb, $ligand_file, $work_dir) = @_;
    my $work_out = "$work_dir/out";

    #
    # Determine if our protein is large enough to require a smaller batch_size
    #
    my $residues;
    my $ok = IPC::Run::run(["count-pdb-residues", $pdb->{local_path}], ">", \$residues);
    chomp $residues;
    if (!$ok)
    {
	die "Could not determine residue coutn for $pdb->{local_path}";
    }
    print STDERR "PDB has $residues residues\n";
    if ($residues > 1024 && !defined($self->params->{batch_size}))
    {
	$self->params->{batch_size} = 5;
	print STDERR "Setting batch_size to " . $self->params->{batch_size} . "\n";
    }
    
    my @batch_size;
    if ($self->params->{batch_size} =~ /\d/)
    {
	@batch_size = ('--batch-size', $self->params->{batch_size});
    }
    my @cmd = ('run_local_docking',
	       @batch_size,
		'--drug-dbs', $ligand_file,
		'--name', 'diffdock_1_1',
		'--receptor-pdb', $pdb->{local_path},
	       $work_out);
    print "RUN @cmd\n";
    #
    # Quasi debugging; don't run if the output is already there.
    if (! -d $work_out)
    {
	my $ok = IPC::Run::run(\@cmd);
	
	print "ok=$ok\n";
	if (!$ok)
	{
	    die "Docking run failed\n";
	}
    }

    #
    # Copy output from work dir to the output directory.
    # We preserve ligand output directory from the docking run
    # and just move it into the output directory.
    #
    # While we're there, we load the result.csv and store the contents
    # for later reporting.
    #

    my $out = $self->output_dir . "/$pdb->{pdb_id}";
    make_path($out);

    for my $ligand (@{$self->ligand_name})
    {
	print "Copy $work_out/$ligand to $out\n";
	if (!dircopy("$work_out/$ligand", "$out/$ligand"))
	{
	    die "Failed to copy $work_out/$ligand to $out";
	}
	#
	# Also copy the diffdock logfile.
	#
	copy("$work_out/diffdock_log", "$out/diffdock_log.txt");
	
	my $result_data = csv(in => "$work_out/$ligand/result.csv", headers => 'auto', sep_char => "\t");
	$_->{output_folder} = "$pdb->{pdb_id}/$ligand" foreach @$result_data;
	$self->{result_data}->{$pdb->{pdb_id}}->{$ligand} = $result_data;
    }
}

sub load_ligand_smiles
{
    my($self, $smiles_list) = @_;

    #
    # We just store this in the smiles_list member.
    #
    $self->{smiles_list} = $smiles_list;
    my $staging_dir = $self->staging_dir;
    
    my $file = $self->staging_dir . "/raw_ligands.smi";
    my $validated_ligands_file = $self->staging_dir . "/ligands.smi";
    open(F, ">", $file) or die "Cannot write $file: $!";
    my $row = 0;
    my @new;
    for my $elt_in (@$smiles_list)
    {
	my $elt;
	my $id;
	if (ref($elt_in) eq 'ARRAY')
	{
	    next if @$elt_in == 0;
	    if (@$elt_in == 1)
	    {
		$elt = $elt_in->[0];
	    }
	    else
	    {
		$id = $elt_in->[0];
		$elt = $elt_in->[1];
	    }
	}
	else
	{
	    $elt = $elt_in;
	}
	$id //= sprintf("ligand-%04d", $row + 1);
	
	$self->{ligand_map}->{$id} = $row;
	$self->{ligand_name}[$row] = $id;
	$self->{ligand_info}[$row] = { id => $id, idx => $row, smiles => $elt };
	print F join("\t", $id, $elt), "\n";
	push(@new, $elt);
	$row++;
    }
    $self->{smiles_list} = \@new;

    close(F);

    #RUN the step to check the inputs
    my @cmd = (
        "check_input_smile_strings",
        "$staging_dir",
        "$file",
        );
    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Input clean up command failed $?: @cmd";
    }

    return $validated_ligands_file;
}


sub load_ligand_library {
    my ($self, $lib_name) = @_;
    #
    # Given a ligand library write ID, Name, and SMILE string to info.txt with names
    # And pass only ID and SMILE string to match the other inputs
    #
    open(my $fh, '<', $lib_name) or die "Could not open file '$lib_name' $!";
    my $file_contents = do { local $/; <$fh> };
    close($fh);

    my $dat = ($file_contents);
    my $staging_dir = $self->staging_dir;

    # Open info.txt for writing in the staging directory
    open(my $info_fh, '>', "$staging_dir/info.txt") or die "Could not open file '$staging_dir/info.txt' $!";

    open(IN, "<", \$dat) or die "Cannot string-open results: $!";
    my @dat;
    while (<IN>) {
        chomp;
        s/^\s*//; # Remove leading whitespace
        next if $_ eq '';
        my @cols = split(/\s+/);
        # Check if there are exactly three columns
        if (scalar @cols == 3) {
            # Push only ID and SMILES string to @dat
            push(@dat, [$cols[0], $cols[2]]); # Only pass the ID and SMILES to @dat

            # Write all three columns to info.txt
            print $info_fh join(' ', @cols) . "\n"; 
        }
    }
    close($info_fh);
    my $res = $self->load_ligand_smiles(\@dat);
    return $res;
}

sub load_ligand_ws_file
{
    my($self, $ws_file) = @_;

    my $dat = $self->app->workspace->download_file_to_string($ws_file);

    open(IN, "<", \$dat) or die "Cannot string-open results: $!";
    my @dat;
    while (<IN>)
    {
	chomp;
	s/^\s*//;
	next if $_ eq '';
	my @cols = split(/\s+/);
	push(@dat, [@cols]);
    }

    my $res = $self->load_ligand_smiles(\@dat);
    return $res;
}

sub stage_pdb
{
    my($self, $pdb_list) = @_;

    my $qry = join(",", @$pdb_list);
    my @res = $self->api->query('protein_structure', ['in', 'pdb_id', '(' . $qry . ')'], ['select', 'pdb_id,gene,product,method,patric_id,title,file_path']);

    my %res = map { ($_->{pdb_id} => $_ ) } @res;

    my @out;
    for my $pdb (@$pdb_list)
    {
	my $ent = $res{$pdb};
	if (!$ent)
	{
	    die "No pdb found for $pdb\n";
	}
	my $path = "/vol/structure/$ent->{file_path}";
	if (! -s $path)
	{
	    die "Cannot find $path\n";
	}
	my $local = $self->staging_dir . "/" . basename($path);
	copy($path, $local);
	$ent->{local_path} = $local;
	$self->{pdb_info}->{$pdb} = $ent;
	push(@out, $ent);
    }
    return @out;
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $mem = '128G';
    
    my $time = 60 * 60 * 10;
    my $pf = {
	cpu => 8,
	memory => $mem,
	runtime => $time,
	policy_data => { gpu_count => 1, partition => 'gpu' },
    };
    return $pf;
}

sub save_output_files
{
    my($self, $app, $output) = @_;
    my %suffix_map = (
		      bai => 'bai',
		      bam => 'bam',
		      #
		      # The result.csv generated by the app is really tab separated.
		      #
		      csv => 'tsv',
		      depths => 'txt',
		      err => 'txt',
		      fasta => "contigs",
		      html => 'html',
		      out => 'txt',
		      txt => 'txt',
		      png => 'png',
		      pdb => 'pdb',
		      tsv => 'tsv',
		      txt => 'txt',);

    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(D, $output))
    {
	while (my $p = readdir(D))
	{
	    next if ($p =~ /^\./);
	    my @cmd = ("p3-cp", "--overwrite", "--recursive", @suffix_map, "$output/$p", "ws:" . $app->result_folder);
	    print STDERR "saving files to workspace... @cmd\n";
	    my $ok = IPC::Run::run(\@cmd);
	    if (!$ok)
	    {
		warn "Error $? copying output with @cmd\n";
	    }
	}
    closedir(D);
    }
}

1;
