function deleteGenomeConfirmation(user,genome,key){
	console.log("$ user      = '"+user+"'");
	console.log("$ project   = '"+genome+"'");
	console.log("$ key       = '"+key+"'");
	panel_iframe             = document.getElementById('panel_genome_iframe');
	dom_object               = panel_iframe.contentDocument.getElementById('g_delete_'+key);

	dom_object.innerHTML     = "<b><font color=\"red\">[Are you sure?]</font><button type='button' onclick='parent.deleteGenome_yes(\""+user+"\",\""+genome+"\",\""+key+"\")'>Yes, delete.</button>";
    dom_object.innerHTML    += "<button type='button' onclick='parent.deleteGenome_no(\""+user+"\",\""+genome+"\",\""+key+"\")'>No, cancel</button></b>";

	dom_button               = panel_iframe.contentDocument.getElementById('genome_delete_'+key);
	dom_button.style.display = 'none';
}

function deleteGenome_yes(user,genome,key){
	$.ajax({
		url : 'genome.delete_server.php',
		type : 'post',
		data : {
			user: user,
			genome: genome
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				// reload entire page - in order to ensure the update of the quota calculation
				window.top.location.reload();
			}
		}
	});
}


function deleteGenome_no(user,genome,key){
	panel_iframe             = document.getElementById('panel_genome_iframe');
	dom_object               = panel_iframe.contentDocument.getElementById('g_delete_'+key)
	dom_object.innerHTML     = "";

	dom_button               = panel_iframe.contentDocument.getElementById('genome_delete_'+key);
	dom_button.style.display = 'inline';
}
