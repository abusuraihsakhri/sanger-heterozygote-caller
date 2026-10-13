"use strict";
(function () {
    const $ = id => document.getElementById(id);
    let trace = null;
    let results = [];
    let currentName = "";
    let isDemo = false;

    function setStatus(message, error) {
        $("status").textContent = message;
        $("status").classList.toggle("error", Boolean(error));
    }
    function settings() {
        return {
            threshold: Number($("threshold").value),
            window: Number($("window").value),
            minQuality: Number($("quality").value)
        };
    }
    function render() {
        if (!trace) return;
        const params = settings();
        results = SangerABIF.callPeaks(trace, params.threshold, params.window, params.minQuality);
        const flagged = results.filter(row => row.flagged);
        $("bases").textContent = String(trace.sequence.length);
        $("reviewed").textContent = String(results.length);
        $("flagged").textContent = String(flagged.length);
        $("sample").textContent = (isDemo ? "Synthetic example: " : "File: ") + currentName;
        const body = $("rows");
        body.replaceChildren();
        if (!flagged.length) {
            const tr = document.createElement("tr");
            const td = document.createElement("td");
            td.colSpan = 7;
            td.textContent = "No candidate peaks at the selected settings.";
            tr.append(td);
            body.append(tr);
        }
        flagged.forEach((row) => {
            const tr = document.createElement("tr");
            const values = [
                row.position, row.calledBase, row.iupac,
                row.primaryBase + ":" + row.primaryHeight + " / " +
                    row.secondaryBase + ":" + row.secondaryHeight,
                row.ratio.toFixed(3), row.quality
            ];
            values.forEach(value => {
                const td = document.createElement("td");
                td.textContent = String(value);
                tr.append(td);
            });
            const td = document.createElement("td");
            const button = document.createElement("button");
            button.className = "secondary";
            button.textContent = "Inspect";
            button.type = "button";
            button.setAttribute("aria-label", "Inspect base " + row.position);
            button.addEventListener("click", () => plot(row));
            td.append(button);
            tr.append(td);
            body.append(tr);
        });
        $("json").disabled = false;
        $("csv").disabled = false;
        if (flagged.length) plot(flagged[0]);
        else {
            const ctx = $("chart").getContext("2d");
            ctx.clearRect(0, 0, $("chart").width, $("chart").height);
            $("chart-note").textContent = "No flagged positions to plot.";
        }
        setStatus("Analyzed " + currentName + ": " + flagged.length + " candidate positions." +
            (isDemo ? " Synthetic example only." : ""));
    }
    function plot(row) {
        const canvas = $("chart");
        const ctx = canvas.getContext("2d");
        const w = canvas.width, h = canvas.height;
        ctx.clearRect(0, 0, w, h);
        const start = Math.max(0, row.traceIndex - 48);
        const end = Math.min(trace.channels.A.length, row.traceIndex + 49);
        let maximum = 1;
        for (const base of ["A", "C", "G", "T"]) {
            for (let x = start; x < end; x++) maximum = Math.max(maximum, trace.channels[base][x]);
        }
        const innerWidth = w - 40, innerHeight = h - 28;
        ctx.strokeStyle = "#506078";
        ctx.beginPath();
        ctx.moveTo(24, 10); ctx.lineTo(24, h - 18); ctx.lineTo(w - 16, h - 18);
        ctx.stroke();
        const colors = {A:"#58d17b", C:"#64b5ff", G:"#e4cf55", T:"#f18684"};
        for (const base of ["A", "C", "G", "T"]) {
            ctx.strokeStyle = colors[base];
            ctx.lineWidth = 2;
            ctx.beginPath();
            const values = trace.channels[base];
            for (let x = start; x < end; x++) {
                const px = 24 + (x - start) / Math.max(1, end - start - 1) * innerWidth;
                const py = h - 18 - Math.max(0, values[x]) / maximum * innerHeight;
                if (x === start) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            }
            ctx.stroke();
        }
        const mx = 24 + (row.traceIndex - start) / Math.max(1, end - start - 1) * innerWidth;
        ctx.strokeStyle = "#f7f7fc";
        ctx.setLineDash([4, 5]);
        ctx.beginPath(); ctx.moveTo(mx, 12); ctx.lineTo(mx, h - 18); ctx.stroke();
        ctx.setLineDash([]);
        ctx.font = "13px sans-serif";
        let labelX = 32;
        for (const base of ["A", "C", "G", "T"]) {
            ctx.fillStyle = colors[base]; ctx.fillText(base, labelX, 20); labelX += 30;
        }
        $("chart-note").textContent =
            "Base " + row.position + " (trace sample " + row.traceIndex + "), " +
            row.primaryBase + "/" + row.secondaryBase + ", ratio " + row.ratio.toFixed(3) +
            ". A green, C blue, G yellow, T red.";
    }
    function save(content, fileName, type) {
        const blob = new Blob([content], {type});
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = fileName;
        document.body.append(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
    }
    function filePrefix() {
        return currentName.replace(/\.(ab1|abi)$/i, "").replace(/[^a-zA-Z0-9._-]/g, "_").slice(0, 80) || "trace";
    }
    function exportJSON() {
        if (!trace) return;
        const payload = {
            sample: currentName, synthetic_example: isDemo, parameters: settings(),
            total_bases: trace.sequence.length, reviewed_bases: results.length,
            flagged_positions: results.filter(row => row.flagged),
            note: "Screening candidates only; not clinically validated."
        };
        save(JSON.stringify(payload, null, 2), filePrefix() + "-candidates.json", "application/json");
    }
    function exportCSV() {
        if (!trace) return;
        const fields = ["position","traceIndex","calledBase","quality","primaryBase","primaryHeight",
            "secondaryBase","secondaryHeight","ratio","iupac","flagged"];
        const lines = [fields.join(",")];
        results.filter(row => row.flagged).forEach(row => {
            lines.push(fields.map(f => String(row[f])).join(","));
        });
        save(lines.join("\r\n") + "\r\n", filePrefix() + "-candidates.csv", "text/csv");
    }
    async function loadFile() {
        const file = $("ab1").files[0];
        if (!file) { setStatus("Choose an AB1 file before analyzing.", true); return; }
        if (file.size > 32 * 1024 * 1024) {
            setStatus("File exceeds the 32 MiB safety limit.", true); return;
        }
        $("analyze").disabled = true;
        setStatus("Reading the selected file locally...");
        try {
            const parsed = SangerABIF.parseAbif(await file.arrayBuffer());
            trace = parsed;
            currentName = file.name;
            isDemo = false;
            render();
        } catch (error) {
            trace = null;
            results = [];
            $("json").disabled = $("csv").disabled = true;
            $("bases").textContent = $("reviewed").textContent = $("flagged").textContent = "—";
            $("rows").replaceChildren();
            setStatus("Unable to analyze file: " + error.message, true);
        } finally {
            $("analyze").disabled = false;
        }
    }
    $("analyze").addEventListener("click", loadFile);
    $("demo").addEventListener("click", () => {
        trace = SangerABIF.demoTrace();
        currentName = "synthetic-demo";
        isDemo = true;
        try { render(); } catch (error) { setStatus(error.message, true); }
    });
    $("json").addEventListener("click", exportJSON);
    $("csv").addEventListener("click", exportCSV);
    for (const key of ["threshold", "window", "quality"]) {
        $(key).addEventListener("change", () => {
            if (!trace) return;
            try { render(); } catch (error) { setStatus(error.message, true); }
        });
    }
})();
