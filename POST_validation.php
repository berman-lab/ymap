<?php
//
// require_once 'POST_validation.php';
// session_destroy();
// header('Location: ../../');
//

function stripHTML_POST($POST_name) {
	// strip out any HTML tags.
	return filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING);
}
function sanitize_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// convert any spaces to underlines.
	$cleanString = str_replace(" ","_", $cleanString);
	// remove whitespace.
	$cleanString = preg_replace("/[\s]+/", "", $cleanString);
	// remove everything but alphanumeric characters, underlines, dashes, and periods.
	$cleanString = preg_replace("/[^\w-_.]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeColor_POST($POST_name) {
	$cleanString = sanitize_POST($POST_name);
	// convert any underlines to spaces.
	$cleanString = str_replace("_"," ", $cleanString);
	return $cleanString;
}
function sanitizeHapmap_POST($POST_name) { //whitelisted characters for hapmap description.
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// convert any spaces to underlines.
	$cleanString = str_replace(" ","_", $cleanString);
	// remove whitespace.
	$cleanString = preg_replace("/[\s]+/", "", $cleanString);
	// remove everything but alphanumeric characters, underlines, dashes, and periods.
	$cleanString = preg_replace("/[^\w-_.\[\]\(\):;]+/", "", $cleanString);
	// convert any underlines to spaces.
	$cleanString = str_replace("_"," ", $cleanString);
	return $cleanString;
}
function sanitizeFile_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// convert any spaces to underlines.
	$cleanString = str_replace(" ","_", $cleanString);
	// remove everything but alphanumeric characters, underlines, periods, and commas (for multiple file strings).
	$cleanString = preg_replace("/[^\w.,]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeFloat_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but numerals and period.
	$cleanString = preg_replace("/[^\d\.]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeInt_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but numerals.
	$cleanString = preg_replace("/[^\d]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeIntChar_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but numerals.
	$cleanString = preg_replace("/[^\d]+/", "", $cleanString);
	// only use first numeral of input.
	$cleanString = $cleanString[0];
	return $cleanString;
}
function sanitizeTabbed_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but word characters and tabs.
	$cleanString = preg_replace("/[^\w\t]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeName_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but letters, spaces, apostraphes, and dashes.
	$cleanString = preg_replace("/[^a-zA-Z-' ]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeEmail_POST($POST_email) {
	if (filter_var($POST_email, FILTER_VALIDATE_EMAIL)) {
		// valid email address.
		return $POST_email;
	} else {
		// invalid email address.
		return "invalid_email_provided";
	}
}

?>
