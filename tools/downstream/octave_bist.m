% Run unmodified upstream BIST blocks using the selected Python environment.
addpath(fullfile(getenv('SYMOX_OCTSYMPY_SOURCE'), 'inst'));
setenv('PYTHON', getenv('SYMOX_AUDIT_PYTHON'));
sympref('quiet', true);
syms x;  % Fail immediately if the Python bridge cannot initialize.
names = {'sym', 'diff', 'int', 'solve'};
rows = {};
failed = 0;
for k = 1:numel(names)
  filename = fullfile(getenv('SYMOX_OCTSYMPY_SOURCE'), 'inst', '@sym', [names{k} '.m']);
  [n, total, xfail, bugs, skipped, runtime_skipped, regressions] = test(filename, 'normal', stdout);
  failures = max(0, total - n - xfail - bugs);
  rows{end+1} = struct('file', names{k}, 'passed', n, 'total', total, ...
    'xfailed', xfail + bugs, 'skipped', skipped + runtime_skipped, ...
    'failed', failures, 'regressions', regressions);
  failed += failures;
endfor
fprintf('SYMOX_BIST_RESULTS=%s\n', jsonencode(rows));
sympref('reset');
exit(failed > 0);
