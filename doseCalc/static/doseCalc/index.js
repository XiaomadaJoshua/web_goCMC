$.validator.addMethod('minStrict', function (value, el, param) {
    return value > param;
}, "Please enter a positive number.");

$( document ).ready(function() {
    $('#CT_upload').click(function(){
		$('#CT_uploader').show();
		$('#material_radio').hide();
	});
	$('#standard_material').click(function(){
		$('#CT_uploader').hide();
		$('#material_radio').show();
	});
	
	$( "#setup_form" ).validate({
		rules: {
			energy: {
				required: true,
				range: [0, 450]
			},
			nhistories:{
				required: true,
				min: 1000
			},
		}
	});
	$('.beamsize').each(function() {
		$(this).rules('add', {
			required: true,
			min: 0
		});
	});
	$('.voxsize').each(function(){
		$(this).rules('add',{
			required: true,
			minStrict: 0
		});
	});
	$('.resolution').each(function() {
		$(this).rules('add', {
			required: true,
			min: 1
		});
	});

	$('#submit').click(function(){
		$('#simulating').show();
	});
});