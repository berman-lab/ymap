<?php
//
// require_once 'POST_validation.php';
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
	// remove everything but alphanumeric characters and underlines.
	$cleanString = preg_replace("/[\s\W]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeFloat_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but numerals and period.
	$cleanString = preg_replace("/[^\d\.]+/", "", $cleanString);
	return $cleanString;
}
function sanitizeIntChar_POST($POST_name) {
	// strip out any html tags.
	$cleanString = trim(filter_input(INPUT_POST, $POST_name, FILTER_SANITIZE_STRING));
	// remove everything but numerals and period.
	$cleanString = preg_replace("/[^\d\]+/", "", $cleanString);
	// only use first numeral of input.
	$cleanString = $cleanString[0];
	return $cleanString;
}


?>
