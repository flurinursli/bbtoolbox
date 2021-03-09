flag = .false.

IF (PRESENT(strict)) flag = strict

IF (PRESENT(v1) .and. PRESENT(v2)) THEN
  IF (flag) THEN
    assert = ANY(par .le. v1) .or. ANY(par .ge. v2)
  ELSE
    assert = ANY(par .lt. v1) .or. ANY(par .gt. v2)
  ENDIF

ELSEIF (PRESENT(v1) .and. .not.PRESENT(v2)) THEN
  IF (flag) THEN
    assert = ANY(par .le. v1)
  ELSE
    assert = ANY(par .lt. v1)
  ENDIF

ELSE
  IF (flag) THEN
    assert = par(1) .ge. par(2)
  ELSE
    assert = par(1) .gt. par(2)
  ENDIF
ENDIF
