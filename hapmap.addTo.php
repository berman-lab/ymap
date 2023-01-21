<?php
	session_start();
	error_reporting(E_ALL);
        require_once 'constants.php';
        ini_set('display_errors', 1);

	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php }
        else {
                $user      = $_SESSION['user'];
                $key       = preg_replace('/\D/', '', $_GET['k']);  //Strip all non-numerical characters from string.
        }
        require_once 'constants.php';
?>
<!DOCTYPE html>
<!--[if lt IE 7]>      <html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>         <html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>         <html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js"> <!--<![endif]-->
<head>
    <meta charset="UTF-8">
    <link rel="stylesheet" href="css/normalize.css">     <!-- Normalizes cross-browser display differences.  --!>
    <link rel="stylesheet" href="css/style.css">         <!-- Normalizes cross-browser styling.              --!>
    <link rel="stylesheet" href="css/defaultTheme.css">  <!-- --!>
	<style type="text/css">
		.tab {
			margin-left:    1cm;
		}
	</style>
	<base target="_parent">
</head>
<BODY class="tab">
<!---------- Main section of User interface ----------!>
	<div class="upload-wrapper">
		<table><tr><td>
		<!--- Button to add haplotype entry to hapmap ---!>
		<button onclick="parent.parent.show_hidden('Hidden_AddToHapmap'); parent.parent.reload_hidden('Hidden_AddToHapmap','hapmap.addTo_window.php?k=<?php echo $key; ?>');">Add haplotype entry...</button>
	</div>
</BODY>
</html>
