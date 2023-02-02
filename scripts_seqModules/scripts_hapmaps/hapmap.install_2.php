<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
	require_once '../../POST_validation.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Load user string from session.
	$user       = $_SESSION['user'];

	// Validate input strings.
	$hapmap          = sanitize_POST("hapmap");
	$genome          = sanitize_POST("genome");
	$project1        = sanitize_POST("project1");
	$project2        = sanitize_POST("project2");
	$referencePloidy = sanitizeFloat_POST("referencePloidy");

	$hapmap_dir = "../../users/".$user."/hapmaps/".$hapmap;

	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if (!(is_dir($genome_dir1) || is_dir($genome_dir2))) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	if ($project1 != '') {
		// Confirm if requested project1 project exists.
		$project1_dir1 = "../../users/".$user."/projects/".$project1;
		$project1_dir2 = "../../users/default/projects/".$project1;
		if (!(is_dir($project1_dir1) || is_dir($project1_dir2))) {
			// Parent project doesn't exist, should never happen: Force logout.
			session_destroy();
			?><script type="text/javascript"> parent.location.reload(); </script><?php
		}
	}

	if ($project2 != '') {
		// Confirm if requested project2 project exists.
		$project2_dir1 = "../../users/".$user."/projects/".$project2;
		$project2_dir2 = "../../users/default/projects/".$project2;
		if (!(is_dir($project2_dir1) || is_dir($project2_dir2))) {
			// Project2 project doesn't exist, should never happen: Force logout.
			session_destroy();
			?><script type="text/javascript"> parent.location.reload(); </script><?php
		}
	}

	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	// Load the number of chromosomes from "chromosome_sizes.txt"
	$file1 = $genome_dir1."/chromosome_sizes.txt";
	$file2 = $genome_dir2."/chromosome_sizes.txt";
	if (file_exists($file1)) {
		$handle     = fopen($file1,'r');
		$fileString = fread($handle, filesize($file1));
	} else {
		$handle     = fopen($file2,'r');
		$fileString = fread($handle, filesize($file2));
	}
	fclose($handle);
	$chrLinesArray  = explode("\n", trim($fileString));
	array_shift($chrLinesArray);
	$chr_count      = count($chrLinesArray);

	$chrNumID       = array();
	$chrNameShort   = array();
	$chrSize        = array();
	$newTags        = array();
	for ($chr=0; $chr<$chr_count; $chr+=1) {
		array_push($chrNumID,"");
		array_push($chrNameShort,"");
		array_push($chrSize,"");
		array_push($newTags,"");
	}
	foreach ($chrLinesArray as $key=>$chrLine) {
		$chrLineParts       = preg_split('/\s+/', $chrLine);
		$oneChrNum          = $chrLineParts[0];
		$oneChrSize         = $chrLineParts[1];
		$oneChrName         = $chrLineParts[2];
		$chrNumID[$key]     = $oneChrNum;
		$chrNameShort[$key] = $oneChrName;
		$chrSize[$key]      = $oneChrSize;
	}

	for ($chr=0; $chr<$chr_count; $chr+=1) {
		$newTags[$chr] = sanitize_POST("newTag_".$chr);
	}
	$newTags[0] = sanitize_POST("newTag");
	print_r($newTags[0]);

	echo "<script type='text/javascript'>\n";
	echo "chr_count = ".$chr_count."\n";
	echo "</script>\n";
?>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
		.tab {margin-left:   1cm;}
		select.list1 option.optionCyan {background-color: #007700;}
	</style>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>[Needs Title]</title>
</HEAD>
<BODY onload="showColors('homolog_a_color','colorBar1'); showColors('homolog_b_color','colorBar2')">
<b>Fill in details for this hapmap entry.</b>
<div class="tab">
	Genome : <?php echo $genome; ?><br>
	<?php
	if ($referencePloidy == 2) {
		echo "Reference dataset : ".$project1."<br>";
	} else {
		echo "Reference dataset 1 : ".$project1."<br>";
	}
	if (file_exists($project1_dir1)) {
	    if (file_exists($project1_dir1."/fig.CNV-LOH-map.2.png")) {
	        $imageUrl = $project1_dir1."/fig.CNV-LOH-map.2.png";
	    } else {
	        $imageUrl = $project1_dir1."/fig.CNV-SNP-map.2.png";
	    }
	    echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
	} else {
		if (file_exists($project1_dir2."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = $project1_dir2."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = $project1_dir2."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
	}
    ?><br>
	<?php
	if ($referencePloidy == 2) {
		echo "Experimental dataset : ".$project2."<br>";
	} else {
		echo "Reference dataset 2 : ".$project2."<br>";
	}
	if (file_exists($project2_dir1)) {
		if (file_exists($project2_dir1."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = $project2_dir1."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = $project2_dir1."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
		$CGD_annotations_url = $project2_dir1."/CGD_annotations.".$project2.".txt";
		if (file_exists($project2_dir1."/CGD_annotations.".$project2.".txt")) {
			if ($referencePloidy == 2) {
				echo "<br><div class='tab'>Examine <button onclick=\"loadExternal('".$CGD_annotations_url."',50,  220);\">GBrowse annotation track</button> to determine precise breakpoints.</div><br>";
			}
		}
	} else {
		if (file_exists($project2_dir2."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = $project2_dir2."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = $project2_dir2."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
		$CGD_annotations_url = $project2_dir2."/CGD_annotations.".$project2.".txt";
		if (file_exists($project2_dir2."/CGD_annotations.".$project2.".txt")) {
			if ($referencePloidy == 2) {
				echo "<button onclick=\"loadExternal('".$CGD_annotations_url."',50,  220);\">GBrowse</button><br>";
			}
		}
	}
	?>
	<br>
</div>
<table><tr><td>
<?php
if ($referencePloidy == 2) {
?>
	<div class="tab">
		<table border="0">
		<tr>
			<th bgcolor=\"#CCFFCC\"><i>Chr</i></th>
			<th bgcolor=\"#CCCCFF\"><i>Homolog</i></th>
			<th bgcolor=\"#CCFFCC\"><i>Start bp</i></th>
			<th bgcolor=\"#CCCCFF\"><i>End bp</i></th>
			<th bgcolor=\"#CCFFCC\"></th>
			<th bgcolor=\"#CCCCFF\"></th>
			<th bgcolor=\"#CCFFCC\"></th>
		</tr>
		<?php
		for ($chr=0; $chr<$chr_count; $chr+=1) {
			$chrID = $chr+1;
			echo "\t\t<tr>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\">{$chrNameShort[$chr]}</td>\n";
			echo "\t\t\t                                         <input type=\"hidden\"   id=\"chr_{$chr}\"     value=\"{$chrNumID[$chr]}\"></input>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><select id=\"homolog_{$chr}\"><option value=\"a\">a</option><option value=\"b\">b</option></select></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><input type=\"text\"     id=\"start_{$chr}\"   value=\"1\"                 size=\"8\"></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><input type=\"text\"     id=\"end_{$chr}\"     value=\"{$chrSize[$chr]}\"  size=\"8\"></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><button onclick=\"updateTag_{$chr}();\">+</button></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><div id=\"tag_{$chr}\">{$newTags[$chr]}</div></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><button onclick=\"clearTag_{$chr}();\">-</button></td>\n";
			echo "\t\t</tr>\n";
		}
		?>
		</table>
		<script type="text/javascript">
		<?php
		for ($chr=0; $chr<$chr_count; $chr+=1) {
			echo "updateTag_{$chr} = function() {\n";
			echo "\tchr_id        = document.getElementById(\"chr_{$chr}\"    ).value;\n";
			echo "\thomolog_id    = document.getElementById(\"homolog_{$chr}\").value;\n";
			echo "\tstart_bp      = document.getElementById(\"start_{$chr}\"  ).value;\n";
			echo "\tend_bp        = document.getElementById(\"end_{$chr}\"    ).value;\n";
			echo "\tnewString     = String(homolog_id)+'['+String(start_bp)+':'+String(end_bp)+'] '\n";
			echo "\tdocument.getElementById(\"tag_{$chr}\").innerHTML += newString;\n";
			echo "}\n";
			echo "clearTag_{$chr} = function() {\n";
	        echo "\tdocument.getElementById(\"tag_{$chr}\").innerHTML = \"\";\n";
	        echo "}\n";
		}
		?>
		function loadExternal(imageUrl) {
			window.open(imageUrl);
		}
		</script>
	</div>
<?php
}
?>
</td><td valign="top">
<?php
if ($referencePloidy == 2) {
?>
    <div class="tab">
        <br>
        If the experimental dataset was analyzed using the reference dataset, any regions of loss of heterozygosity (LOH) would appear in red.<br>
        If the strains are unrelated, it is unlikely this analysis will be helpful.<br><br>
        Enter a description of the homolog identity of the homozygosed regions using the form at left.<br>
        <div class="tab">
            Set the '<i><b>Homolog</b></i>', '<i><b>Start bp</b></i>', and '<i><b>End bp</b></i>' options to describe a LOH region on a chromosome.<br>
            Click the <b>[+]</b> button to accept the LOH description.<br>
            Click the <b>[-]</b> buttons to clear the entry for a chromosome.
        </div>
		<br>
		Clicking the '<b>GBrowse</b>' button will open a custom GBrowse annotation track file for the experimental dataset in a new window.   The
		annotation file can be loaded into GBrowse at <a href="http://Candidagenome.org">Candida Genome Database</a>, for supported genomes.   This
		can be useful in determining the end-coordinates of LOH regions.
    </div>
<?php
}
?>
</td></tr>
<tr><td>
<div class="tab">
	<?php
		$colorsFile = $hapmap_dir."/colors.txt";
		if (file_exists($colorsFile)) {
			$handle       = fopen($colorsFile,'r');
			$colorString1 = trim(fgets($handle));
			$colorString2 = trim(fgets($handle));
			?>
			<br>
			Homolog 'a' color : <input type="text" id="homolog_a_color" name="homolog_a_color" value="<?php echo $colorString1; ?>" readonly>
			<div id="colorBar1"></div>
			<br>
			Homolog 'b' color : <input type="text" id="homolog_b_color" name="homolog_b_color" value="<?php echo $colorString2; ?>" readonly>
			<div id="colorBar2"></div>
			<br><br>
			<?php
		} else {
			?>
			<br>
			Choose homolog 'a' color :
			<select id="homolog_a_color" name="homolog_a_color" onchange="showColors('homolog_a_color','colorBar1')">
				<option value="deep pink">deep pink</option>
				<option value="magenta">magenta</option>
				<option value="electric indigo">electric indigo</option>
				<option value="blue">blue</option>
				<option value="dodger blue">dodger blue</option>
				<option value="cyan" selected>cyan</option>
				<option value="spring green">spring green</option>
				<option value="green">green</option>
				<option value="chartreuse">chartreuse</option>
				<option value="yellow">yellow</option>
				<option value="dark orange">dark orange</option>
			</select><div id="colorBar1"><div class="tab" style="background-color:rgb(0  ,255,255)"> &nbsp;</div></div>
			<br>
			Choose homolog 'b' color :
			<select id="homolog_b_color" name="homolog_b_color" onchange="showColors('homolog_b_color','colorBar2')">
				<option value="deep pink">deep pink</option>
				<option value="magenta" selected>magenta</option>
				<option value="electric indigo">electric indigo</option>
				<option value="blue">blue</option>
				<option value="dodger blue">dodger blue</option>
				<option value="cyan">cyan</option>
				<option value="spring green">spring green</option>
				<option value="green">green</option>
				<option value="chartreuse">chartreuse</option>
				<option value="yellow">yellow</option>
				<option value="dark orange">dark orange</option>
			</select><div id="colorBar2"><div class="tab" style="background-color:rgb(255,0  ,255)"> &nbsp;</div></div>
			<br><br>
		<?php
		}
	?>
<div>
<script type="text/javascript">
	showColors=function(listToLookAt, targetToChange) {
		var selectedColor = document.getElementById(listToLookAt).value;
		var select        = document.getElementById(targetToChange);
		select.innerHTML  = '';
		if        (selectedColor == "deep pink")       { select.innerHTML = '<div class="tab" style="background-color:rgb(255,0  ,127)"> &nbsp;</div>';
		} else if (selectedColor == "magenta")         { select.innerHTML = '<div class="tab" style="background-color:rgb(255,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "electric indigo") { select.innerHTML = '<div class="tab" style="background-color:rgb(127,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "blue")            { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "dodger blue")     { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,127,255)"> &nbsp;</div>';
		} else if (selectedColor == "cyan")            { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,255)"> &nbsp;</div>';
		} else if (selectedColor == "spring green")    { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,127)"> &nbsp;</div>';
		} else if (selectedColor == "green")           { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "chartreuse")      { select.innerHTML = '<div class="tab" style="background-color:rgb(127,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "yellow")          { select.innerHTML = '<div class="tab" style="background-color:rgb(255,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "dark orange")     { select.innerHTML = '<div class="tab" style="background-color:rgb(255,127,0  )"> &nbsp;</div>';
		} else {                                         select.innerHTML = '<div class="tab" style="background-color:rgb(127,127,127)"> &nbsp;</div>';
		}
	}
</script>
</td><td valign="top">
<?php
	if (file_exists($colorsFile)) {
	} else {
		?>
		<div class="tab">
			<br>
			Select which colors to be used in displaying the homologs for this haplotype map.<br>
			Colors should be chosen to minimize visual confusion.<br>
			<br>
			Red is not available, as it is the color used to highlight LOHs in regions without haplotype definitions.
		</div>
		<?php
	}
?>
</td></tr>
<tr><td>
<div class="tab">
	<button onclick="GatherAndSubmitData();">Process haplotype entry...</button>
	<script type="text/javascript">
		GatherAndSubmitData=function() {
<?php
if ($referencePloidy == 2) {
?>
			var string_allData = "";
			for (var chr=0;chr<chr_count; chr++) {
				var string_start   = String(chr)+"(";
				var string_start   = document.getElementById("chr_"+String(chr)).value+"(";
				var string_inner   = trim(document.getElementById("tag_"+String(chr)).innerHTML);
				var string_end     = ")";
				if (string_inner.length > 0) {
					string_allData += string_start+string_inner+string_end+"; ";
				}
			}
			var string_final = trim(string_allData);
<?php
}
?>
			var colorA = document.getElementById("homolog_a_color").value;
			var colorB = document.getElementById("homolog_b_color").value;

			var autoSubmitForm = document.createElement('form');
			autoSubmitForm.setAttribute('method','post');
			autoSubmitForm.setAttribute('action','hapmap.install_3.php');
<?php
if ($referencePloidy == 2) {
?>
			var input1 = document.createElement('input');
			    input1.setAttribute('type','hidden');
			    input1.setAttribute('name','hapmap_description');
			    input1.setAttribute('value',string_final);
			    autoSubmitForm.appendChild(input1);
<?php
}
?>
			var input2 = document.createElement('input');
			    input2.setAttribute('type','hidden');
			    input2.setAttribute('name','hapmap');
			    input2.setAttribute('value','<?php echo $hapmap; ?>');
			    autoSubmitForm.appendChild(input2);
			var input3 = document.createElement('input');
			    input3.setAttribute('type','hidden');
			    input3.setAttribute('name','genome');
			    input3.setAttribute('value','<?php echo $genome; ?>');
			    autoSubmitForm.appendChild(input3);
			var input4 = document.createElement('input');
			    input4.setAttribute('type','hidden');
			    input4.setAttribute('name','referencePloidy');
			    input4.setAttribute('value','<?php echo $referencePloidy; ?>');
			    autoSubmitForm.appendChild(input4);
			var input5 = document.createElement('input');
			    input5.setAttribute('type','hidden');
			    input5.setAttribute('name','project1');
			    input5.setAttribute('value','<?php echo $project1; ?>');
			    autoSubmitForm.appendChild(input5);
			var input6 = document.createElement('input');
			    input6.setAttribute('type','hidden');
			    input6.setAttribute('name','project2');
			    input6.setAttribute('value','<?php echo $project2; ?>');
			    autoSubmitForm.appendChild(input6);
			var input7 = document.createElement('input');
			    input7.setAttribute('type','hidden');
			    input7.setAttribute('name','homolog_a_color');
			    input7.setAttribute('value',colorA);
			    autoSubmitForm.appendChild(input7);
			var input8 = document.createElement('input');
			    input8.setAttribute('type','hidden');
			    input8.setAttribute('name','homolog_b_color');
			    input8.setAttribute('value',colorB);
			    autoSubmitForm.appendChild(input8);
			document.body.appendChild(autoSubmitForm);
			autoSubmitForm.submit();
		}
		function trim(str) {
			return str.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
		}
	</script>
</div>
</td><td valign="top">
	<div class="tab">
		<?php
		if ($referencePloidy == 2) {
			echo "Once you've described any LOH regions for this experimental dataset, click the ";
			echo "'<i><b>Process haplotype entry...</b></i>' button to process this entry into the haplotype map.";
		} else {
			echo "Once you've chosen the colors to use for this haplotype map, click the ";
			echo "'<i><b>Process haplotype entry...</b></i>' button to process these datasets into the haplotype map.";
		}
		?>
	</div>
</td><tr></table>
</BODY>
</HTML>
