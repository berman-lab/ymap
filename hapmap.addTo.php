<?php
	session_start();
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
		<table><tr>
		<!--- Button to add haplotype entry to hapmap ---!>
		<button onclick="AddHaplotypeEntry();">Add haplotype entry...</button>

		<script type="text/javascript">
			AddHaplotypeEntry=function() {
			var autoSubmitForm = document.createElement('form');
			autoSubmitForm.setAttribute('method','post');
			autoSubmitForm.setAttribute('action','hapmap.addTo_1.php');
			var input1 = document.createElement('input');
				input1.setAttribute('type','hidden');
				input1.setAttribute('name','user');
				input1.setAttribute('value',user);
				autoSubmitForm.appendChild(input1);
			var input2 = document.createElement('input');
				input2.setAttribute('type','hidden');
				input2.setAttribute('name','hapmap');
				input2.setAttribute('value',hapmap);
				autoSubmitForm.appendChild(input2);
			var input3 = document.createElement('input');
				input3.setAttribute('type','hidden');
				input3.setAttribute('name','key');
				input3.setAttribute('value',key);
				autoSubmitForm.appendChild(input3);
			autoSubmitForm.submit();
			}
		</script>
	</div>
</BODY>
</html>
