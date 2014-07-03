<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">

<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Register New User</title>
	</head>
	<body>
		<form action="php/register_server.php" method="post">
			<label for="user">Username: </label>                                        <input type="text"     id="user"                      name="user"><br>
			<label for="pw1">Password: </label>                                         <input type="password" id="pwOrig"                    name="pwOrig"><br>
			<label for="pw2">Retype Password: </label>                                  <input type="password" id="pwCopy"                    name="pwCopy"><br>
			<b>Optional:</b><br>
			<label for="primaryInvestigator_name">Primary Investigator Name: </label>   <input type="text"     id="primaryInvestigator_name"  name="primaryInvestigator_name"><br>
			<label for="primaryInvestigator_email">Primary Investigator Email: </label> <input type="text"     id="primaryInvestigator_email" name="primaryInvestigator_email"><br>
			<label for="researchInstitution">Research institution: </label>             <input type="text"     id="researchInstitution"       name="researchInstitution"><br>
			<label for="secondary_name">Secondary Name: </label>                        <input type="text"     id="secondary_name"            name="secondary_name"><br>
			<label for="secondary_email">Secondary Email: </label>                      <input type="text"     id="secondary_email"           name="secondary_email"><br>

			<button type="submit">Register</button>
		</form>
		<button type="button" onclick="window.location.href='index.php'">Back to Home</button>
	</body>
</html>
