$(function() {

	// Initialize pager
	var pager = $('#pager');
	pager.pager('init', [ 'Description', 'idock Job', 'Submitted on', 'Status', 'Progress', 'Result' ], function(job) {
		var status, progress, result = '<a href="iview/?' + job._id + '"><img src="/iview/logo.png" alt="iview"></a>';
		if (!job.scheduled) {
			status = 'Queued for execution';
			progress = 0;
		} else if (!job.done) {
			status = 'Execution in progress';
			progress = 0//job.generation / job.num_generations;
		} else {
			status = 'Done on ' + $.format.date(new Date(job.done), 'yyyy/MM/dd HH:mm:ss');
			progress = 1;
			result += '<a href="jobs/' + job._id + '/log.csv.gz"><img src="/excel.png" alt="log.csv.gz"></a><a href="jobs/' + job._id + '/ligands.pdbqt.gz"><img src="/molecule.png" alt="ligands.pdbqt.gz"></a>';
		}
		return [
			job.description,
			'<a href="/idock/iview/?' + job.idock_id + '"><img src="/iview/logo.png" alt="iview"></a>',
			$.format.date(new Date(job.submitted), 'yyyy/MM/dd HH:mm:ss'),
			status,
			(100 * progress).toFixed(5) + '%',
			result
		];
	});

	// Refresh the table of jobs and its pager every second
	var jobs = [], skip = 0;
	var tick = function() {
		$.get('jobs', { skip: skip, count: jobs.length }, function(res) {
			if (res.length) {
				for (var i = skip; i < jobs.length; ++i) {
					var job = res[i - skip];
					jobs[i].done = job.done;
				}
				pager.pager('refresh', skip, jobs.length, 3, 6, false);
				if (res.length > jobs.length - skip) {
					var len = jobs.length;
					jobs = jobs.concat(res.slice(jobs.length - skip));
					pager.pager('source', jobs);
					pager.pager('refresh', len, jobs.length, 0, 6, true);
				}
				for (; skip < jobs.length && jobs[skip].done; ++skip);
			}
			setTimeout(tick, 1000);
		});
	};
	tick();

	// Initialize sliders
	$('#mms').slider({
		range: true,
		min: 0,
		max: 1000,
		values: [ 300, 500 ]
	});
	$('#nrb').slider({
		range: true,
		min: 0,
		max: 35,
		values: [ 0, 10 ]
	});
	$('#hbd').slider({
		range: true,
		min: 0,
		max: 20,
		values: [ 0, 5 ]
	});
	$('#hba').slider({
		range: true,
		min: 0,
		max: 18,
		values: [ 0, 10 ]
	});
	$('#nha').slider({
		range: true,
		min: 1,
		max: 100,
		values: [ 20, 70 ]
	});
	$('#lgp').slider({
		range: true,
		min: -6,
		max: 12,
		values: [ -0.4, 5.6 ]
	});
	$('#psa').slider({
		range: true,
		min: 0,
		max: 300,
		values: [ 0, 140 ]
	});
	$('#mrf').slider({
		range: true,
		min: 0,
		max: 300,
		values: [ 40, 130 ]
	});
	['mms', 'nrb', 'hbd', 'hba', 'nha', 'lgp', 'psa', 'mrf'].forEach(function(key) {
		var values = $('#' + key).slider('option', 'values');
		$('#' + key + '_lb').text(values[0]);
		$('#' + key + '_ub').text(values[1]);
	});
	$('.slider').slider({
		slide: function(event, ui) {
			$('#' + this.id + '_lb').text(ui.values[0]);
			$('#' + this.id + '_ub').text(ui.values[1]);
		},
	});

	// Initialize tooltips
	$('.form-group a').tooltip();

	// Process submission
	var submit = $('#submit');
	submit.click(function() {
		// Hide tooltips
		$('.form-group a').tooltip('hide');
		// Disable the submit button for a while
		submit.prop('disabled', true);
		// Post a new job with server side validation
		$.post('jobs', {
			idock_id: idockJobId,
			description: $('#description').val(),
			email: $('#email').val(),
			mms_lb: $('#mms_lb').text(),
			mms_ub: $('#mms_ub').text(),
			nrb_lb: $('#nrb_lb').text(),
			nrb_ub: $('#nrb_ub').text(),
			hbd_lb: $('#hbd_lb').text(),
			hbd_ub: $('#hbd_ub').text(),
			hba_lb: $('#hba_lb').text(),
			hba_ub: $('#hba_ub').text(),
			nha_lb: $('#nha_lb').text(),
			nha_ub: $('#nha_ub').text(),
			lgp_lb: $('#lgp_lb').text(),
			lgp_ub: $('#lgp_ub').text(),
			psa_lb: $('#psa_lb').text(),
			psa_ub: $('#psa_ub').text(),
			mrf_lb: $('#mrf_lb').text(),
			mrf_ub: $('#mrf_ub').text(),
		}, function(res) {
			var keys = Object.keys(res);
			// If server side validation fails, show the tooltips
			if (keys.length) {
				keys.forEach(function(key) {
					$('#' + key + '_label').tooltip('show');
				});
			} else {
				$('html, body').animate({ scrollTop: pager.offset().top });
			}
		}, 'json').always(function() {
			submit.prop('disabled', false);
		});
	});

	// Apply accordion to tutorials
	$('.ui-accordion').accordion({
		collapsible: true,
		active: false,
		heightStyle: 'content',
		activate: function(event, ui) {
			$('img', this).trigger('expand');
		}
	});
	$('.ui-accordion img').lazyload({
		event: 'expand',
		effect: 'fadeIn',
	});

	var idockJobId = location.search.substr(1);
	var v = new validator({ id: idockJobId });
	if (v
		.field('id').message('must be a valid object id').objectid().copy()
		.failed()) {
		submit.prop('disabled', true);
		return;
	};
	$.get('/idock/job', {
		id: idockJobId,
	}, function(res) {
		if (res) {
			Object.keys(res).forEach(function(key) {
				$('#' + key).val(res[key]);
			});
		} else {
			submit.prop('disabled', true);
		}
	}, 'json');
});
