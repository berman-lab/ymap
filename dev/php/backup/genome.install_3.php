<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<HTML>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
		.upload {
			width:          675px;      // 675px;
			border:         0;
			height:         40px;   // 40px;
			vertical-align: middle;
			align:          left;
			margin:         0px;
			overflow:       hidden;
		}
		html, body {
			margin:         0px;
			border:         0;
			overflow:       hidden;
		}
	</style>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1">
<title>Install genome into pipeline.</title>
</HEAD>
<?php
    require_once 'constants.php';

	$user             = $_SESSION['user'];
	$key              = filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING);
	$genome           = $_SESSION['genome_'.$key];
//	$chr_count        = $_SESSION['chr_count_'.$key];
//	$chr_used_count   = $_SESSION['chr_used_count_'.$key];
//	$chr_lengths      = $_SESSION['chr_lengths_'.$key];
//	$chr_names        = $_SESSION['chr_names_'.$key];
//	$chr_draws        = $_SESSION['chr_draws_'.$key];
//	$chr_shortNames   = $_SESSION['chr_shortNames_'.$key];
//	$chr_cenStarts    = $_SESSION['chr_cenStarts_'.$key];
//	$chr_cenEnds      = $_SESSION['chr_cenEnds_'.$key];
	$rDNA_chr         = $_SESSION['rDNA_chr_'.$key];
	$rDNA_start       = $_SESSION['rDNA_start_'.$key];
	$rDNA_end         = $_SESSION['rDNA_end_'.$key];
//	$ploidyDefault    = $_SESSION['ploidyDefault_'.$key];
	$annotation_count = $_SESSION['annotation_count_'.$key];

	$annotation_chrs       = array();
	$annotation_shapes     = array();
	$annotation_starts     = array();
	$annotation_ends       = array();
	$annotation_names      = array();
	$annotation_fillColors = array();
	$annotation_edgeColors = array();
	$annotation_sizes      = array();

// Open 'process_log.txt' file.
    $logOutputName = $directory."users/".$user."/genomes/".$genome."/process_log.txt";
    $logOutput     = fopen($logOutputName, 'a');
    fwrite($logOutput, "Running 'php/genome.install_3.php'.\n");

// process POST data.
	fwrite($logOutput, "\tProcessing POST data containing annotation specifications.\n");
    for ($annotation=0; $annotation<$annotation_count; $annotation += 1) {
		$annotation_chr       = $_POST['annotation_chr_'.$annotation];
		$annotation_shape     = $_POST['annotation_shape_'.$annotation];
		$annotation_start     = $_POST['annotation_start_'.$annotation];
		$annotation_end       = $_POST['annotation_end_'.$annotation];
		$annotation_name      = $_POST['annotation_name_'.$annotation];
		$annotation_fillColor = $_POST['annotation_fillColor_'.$annotation];
		$annotation_edgeColor = $_POST['annotation_edgeColor_'.$annotation];
		$annotation_size      = $_POST['annotation_size_'.$annotation];

		$annotation_chrs[$annotation]       = $annotation_chr;
		$annotation_shapes[$annotation]     = $annotation_shape;
		$annotation_starts[$annotation]     = $annotation_start;
		$annotation_ends[$annotation]       = $annotation_end;
		$annotation_names[$annotation]      = $annotation_name;
		$annotation_fillColors[$annotation] = $annotation_fillColor;
		$annotation_edgeColors[$annotation] = $annotation_edgeColor;
		$annotation_sizes[$annotation]      = $annotation_size;
    }

	// Generate 'annotations.txt' :
	fwrite($logOutput, "\tGenerating 'annotations.txt' file.\n");
    $outputName       = $directory."users/".$user."/genomes/".$genome."/annotations.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
	    $output       = fopen($outputName, 'w');
	    fwrite($output, "# Chr\tType\tstart(bp)\tend(bp)\tName\tFillColor\tEdgeColor\tSize\n");
	    if ($rDNA_chr != "null") {
	        fwrite($output, $rDNA_chr."\tdot\t".$rDNA_start."\t".$rDNA_end."\trDNA\tb\tk\t5\n");
	    }
		for ($annotation=0; $annotation<$annotation_count; $annotation += 1) {
			fwrite($output, $annotation_chrs[$annotation]."\t".$annotation_shapes[$annotation]."\t".$annotation_starts[$annotation]."\t".$annotation_ends[$annotation]."\t".$annotation_names[$annotation]."\t".$annotation_fillColors[$annotation]."\t".$annotation_edgeColors[$annotation]."\t".$annotation_sizes[$annotation]."\n");
		}
	}
	fclose($output);

// Debugging output of all variables.
//	print_r($GLOBALS);
?>
<BODY onload = "parent.resize_iframe('<?php echo $key; ?>', 100); document.processing_form.submit();">
	<form name="processing_form" id="processing_form" action="genome.install_4.php" method="post">
	<input type="submit" id="form_submit" name="form_submit" value="Continue to processing..." style="visibility:hidden" disabled>
	<input type="hidden" id="key" name="key" value="<?php echo $key; ?>">
	</form>
	<font color="red">[Installation in process.]</font><br>
</BODY>
</HTML>
<?php
	fwrite($logOutput, "\t'php/genome_install_3.php' has completed.\n");
	fclose($logOutput);
?>
