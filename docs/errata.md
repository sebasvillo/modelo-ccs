# Errata: numerical changes versus the thesis

Every change that moves a number reported by the thesis or by the legacy notebook gets one entry here,
plus an updated test. Golden masters are never rebaselined silently, and numbers are never tuned back
toward the thesis values. Audit references point to [`audit-2026-09.md`](audit-2026-09.md).

| id | equation / module | old value | new value | cause | test | PR |
|----|-------------------|-----------|-----------|-------|------|----|
| E-000 | *(template)* `absorber.py`, Eq. 2.xx | H = … m | H = … m | audit §… : short reason | `tests/…::test_…` | #… |

<!--
Template for a new entry (copy the row, increment the id):
| E-001 | <equation id or module> | <value + units, thesis/legacy> | <value + units> | <cause, audit §> | <tests/path::test_name> | <PR link> |
-->
