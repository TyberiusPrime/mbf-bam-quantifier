---
bookCollapseSection: true
weight: 15
---

# Filter section

We can filter reads using the filters listebd below.

Note the double `[[' in ``[filter]]``, since you're defining a list of filters.

Filters are interpreted eagerly - they're applied as early as possible to save on computation.

All filters have an `action` property that controls whether we restrict
to just reads matching the filter (`keep_only`) or remove those matching the filter( `remove`).

{{<mynav>}}
